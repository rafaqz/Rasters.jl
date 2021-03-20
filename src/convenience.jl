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
- `metadata`: Metadata as a [`StackMetadata`](@ref) object.
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

If the source can't be save as a stack-like object, individual array layers will be saved.
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
    series(dirpath::AbstractString, dims; kw...) => AbstractGeoSeries

Load a vector of filepaths as a `AbstractGeoSeries`. `kw` are passed to the constructor.

`dims` Dimensions can hold an index matching the files in the directory,
or a function to convert the filename strings to index values.
"""
function series end

function series(dirnames; 
    child=geoarray, child_kwargs=nothing, window=(), kw...
)
    filepaths = readdir(path)
    if all(x -> splitext(x)[2], filenames) == splitext(first(filenames))[2]
        # All the same kind of file. We don't need to load them up front.
        DiskStack(filenames; childtype=_get_constructor(geoarray, first(filenames)), kw...)
    else
        # These files are different extensions, just load them all
        # as separate `AbstarctGeoArray` (which has some up front cost from
        # reading the dimensions). They are probably still disk-backed for the actual array.
        arrays = map(filenames) do dn
            _constructor(geoarray, fn)(fn; window=window, child_kwargs...)
        end
        GeoStack(arrays; kw...)
    end
end

# Support methods

function _constructor(method::Function, filename; throw=true)
    _, extension = splitext(filename)
    return if extension in (".grd", ".gri")
        if method === geoarray
            GRDarray
        elseif method === stack
            throw ? _no_stack_error(extension) : nothing
        end
    elseif extension == ".nc"
        _check_imported(:NCDatasets, :NCDarray, extension)
        if method === geoarray
            NCDarray
        else
            NCDstack
        end
    elseif extension == ".h5"
        # In future we may need to examine the file and check if
        # it's a SMAP file or something else that uses .h5
        _check_imported(:HDF5, :SMAPstack, extension)
        if method === geoarray
            throw ? _no_gearray_error(extension) : nothing
        elseif method === stack
            SMAPstack
        end
    else # GDAL handles too many extensions to list, so just try it and see if it works
        _check_imported(:ArchGDAL, :GDALarray, extension)
        if method === geoarray
            GDALarray
        elseif method === stack
            throw ? _no_stack_error(extension) : nothing
        end
    end
end

_no_stack_error(ext) =
    error("$ext files do not have named layers. Use `geoarray(filename)` to load, or `stack((key1=fn1, key2=fn2, ...)`")

_no_gearray_error(ext) =
    error("$ext files not have a single-layer implementation. Use `stack(filename)` to load")

_check_imported(modulename, type, extension)  =
    type in names(GeoData) || error("Run `import $modulename` to enable loading $extension files.")

function series(dirpath::AbstractString, dims=Dim{:series}(); ext=nothing, child=geoarray, kw...)
    filepaths = filter_ext(dirpath, ext)
    series(filepaths, dims; child=child, kw...)
end
function series(filepaths::AbstractVector{<:AbstractString}, dims=Dim{:series}(); child=geoarray, kw...)
    childtype = _constructor(child, first(filepaths))
    GeoSeries(filepaths, dims; childtype=childtype, kw...)
end
