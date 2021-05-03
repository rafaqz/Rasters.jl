"""
    geoarray(filename; kw...) => AbstractGeoArray

Load a file path as `AbstractGeoArray`. 

# Keywords

Passed to the constructor for the file type, and commmonly include:

- `crs`: crs to use instead of the detected crs.
- `mappedcrs`: CRS format like `EPSG(4326)` used in `Selectors` like `Between` and `At`, and
    for plotting. Can be any CRS `GeoFormat` from GeoFormatTypes.jl, like `WellKnownText`.
- `name`: `Symbol` name for the array.
- `dims`: `Tuple` of `Dimension`s for the array. Detected automatically, but can be passed in.
- `refdims`: `Tuple of` position `Dimension`s the array was sliced from.
- `missingval`: Value reprsenting missing values. Detected automatically when possible, but
    can be passed it.
- `metadata`: `Metadata` object for the array. Detected automatically but can be passed in.
"""
function geoarray end


geoarray(filename::AbstractString; kw...) = _constructor(geoarray, filename)(filename; kw...)

"""
    Base.write(filename::AbstractString, A::AbstractGeoArray; kw...)

Write any [`AbstractGeoArray`](@ref) to file, guessing the backend from the file extension.

Keyword arguments are passed to the `write` method for the backend.
"""
function Base.write(
    filename::AbstractString, A::AbstractGeoArray; kw...
)
    write(filename, _constructor(geoarray, filename), A; kw...)
end

"""
    stack(filename; kw...) => AbstractGeoStack
    stack(filenames::NamedTuple; kw...) => AbstractGeoStack

Load a file path or a `NamedTuple` of paths as an `AbstractGeoStack`.

# Arguments

- `filename`: A `NamedTuple` of stack keys and `String` filenames.

# Keywords

Passed to the constructor for the file type, and commmonly include:

- `window`: A `Tuple` of `Dimension`/`Selector`/indices that will be applied to the
    contained arrays when they are accessed.
- `metadata`: `Metadata` as object.
- `child_kwargs`: A `NamedTuple` of keyword arguments to pass to the `childtype` constructor.
- `refdims`: `Tuple` of  position `Dimension` the array was sliced from.

```julia
files = (:temp="temp.tif", :pressure="pressure.tif", :relhum="relhum.tif")
stack = stack(files; child_kwargs=(mappedcrs=EPSG(4326),))
stack[:relhum][Lat(Contains(-37), Lon(Contains(144))
```
"""
function stack end

stack(filename::AbstractString; kw...) = _constructor(stack, filename)(filename; kw...)
stack(filenames::Union{Tuple,AbstractArray}; keys, kw...) = 
    stack(NamedTuple{map(Symbol, keys)}(filenames); kw...) 
function stack(filenames::NamedTuple; child_kwargs=(), kw...)
    if all(x -> splitext(x)[2] == splitext(first(filenames))[2], filenames)
        # All the same kind of file. We don't need to load them up front.
        DiskStack(filenames; childtype=_constructor(geoarray, first(filenames)), kw...)
    else
        # These files are different extensions, just load them all
        # as separate `AbstarctGeoArray` (which has some up front cost from
        # reading the dimensions). They are probably still disk-backed for the
        # actual array.
        arrays = map(filenames) do fn
            _constructor(geoarray, fn)(fn; child_kwargs...)
        end
        GeoStack(arrays; kw...)
    end
end

"""
    Base.write(filename::AbstractString, T::Type{<:AbstractGeoArray}, s::AbstractGeoStack)

Write any [`AbstractGeoStack`](@ref) to file, guessing the backend from the file extension.

Keyword arguments are passed to the `write` method for the backend.

If the source can't be saved as a stack-like object, individual array layers will be saved.
"""
function Base.write(filename::AbstractString, s::AbstractGeoStack; kw...)
    base, ext = splitext(filename)
    T = _constructor(stack, filename; throw=false)
    if T isa Nothing
        # Save as individual `geoarray`
        T = _constructor(geoarray, filename)
        for key in keys(s)
            fn = joinpath(string(base, "_", key, ext))
            write(fn, T, s[key])
        end
    else
        write(filename, _constructor(stack, filename), s; kw...)
    end
end

"""
    series(dirpath::AbstractString, dims; ext, child=geoarray, kw...) => AbstractGeoSeries
    series(filenames::AbstractVector{String}, dims; child=geoarray, kw...) => AbstractGeoSeries

Load a whole folder vector of filepaths as a `AbstractGeoSeries`. 
`kw` are passed to the constructor.

# Arguments

- `dims`: A tuple of dimensions, possibly holding an index matching the files in the directory. 
    By default it will be a numbered `Dim{:series}`.

# Keywords

- `child`: function to load the child object - may be `geoarray` or `stack`, 
    `geoarray` by default.
- `ext`: an extension for the files in a directory. If empty, all files are loaded, 
    and should be the same type.
"""
function series end

function series(dirpath::AbstractString, dims=(Dim{:series}(),); ext=nothing, child=geoarray, kw...)
    filepaths = filter_ext(dirpath, ext)
    series(filepaths, dims; child=child, kw...)
end
function series(filepaths::AbstractVector{<:AbstractString}, dims=(Dim{:series}(),); child=geoarray, kw...)
    childtype = _constructor(child, first(filepaths))
    GeoSeries(filepaths, dims; childtype=childtype, kw...)
end

# Support methods

const EXT = (GRD=(".grd", ".gri"), NCD=".nc", SMAP=".h5")

# The the constructor for a geoarray or stack, based on the
# filename extension. GDAL is the fallback for geoarray as it 
# handles so many file types.
function _constructor(method::typeof(geoarray), filename; throw=true)
    _, extension = splitext(filename)
    return if extension in EXT[:GRD]
        GRDarray
    elseif extension == EXT[:NCD]
        _check_imported(:NCDatasets, :NCDarray, extension; throw=throw)
        NCDarray
    elseif extension == EXT[:SMAP]
        # In future we may need to examine the file and check if
        # it's a SMAP file or something else that uses .h5
        throw ? _no_gearray_error(extension) : nothing
    else # GDAL handles too many extensions to list, so just try it and see if it works
        _check_imported(:ArchGDAL, :GDALarray, extension; throw=throw)
        GDALarray
    end
end
function _constructor(method::typeof(stack), filename; throw=true)
    _, extension = splitext(filename)
    return if extension in EXT[:GRD]
        throw ? _no_stack_error(extension) : nothing
    elseif extension == EXT[:NCD]
        _check_imported(:NCDatasets, :NCDarray, extension; throw=throw)
        NCDstack
    elseif extension == EXT[:SMAP]
        _check_imported(:HDF5, :SMAPstack, extension; throw=throw)
        SMAPstack
    else
        throw ? _no_stack_error(extension) : nothing
    end
end

_no_stack_error(ext) =
    error("$ext files do not have named layers. Use `geoarray(filename)` to load, or `stack((key1=fn1, key2=fn2, ...)`")

_no_gearray_error(ext) =
    error("$ext files not have a single-layer implementation. Use `stack(filename)` to load")

function _check_imported(modulename, type, extension; throw=true)
    if type in names(GeoData)
        true
    else
        if throw
            error("Run `import $modulename` to enable loading $extension files.")
        else
            false
        end
    end
end
