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

geoarray(filename::AbstractString; kw...) = GeoArray(filename; kw...)

"""
    Base.write(filename::AbstractString, A::AbstractGeoArray; kw...)

Write any [`AbstractGeoArray`](@ref) to file, guessing the backend from the file extension.

Keyword arguments are passed to the `write` method for the backend.
"""
function Base.write(
    filename::AbstractString, A::AbstractGeoArray; kw...
)
    write(filename, _sourcetype(filename), A; kw...)
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
- `refdims`: `Tuple` of  position `Dimension` the array was sliced from.

```julia
files = (:temp="temp.tif", :pressure="pressure.tif", :relhum="relhum.tif")
stack = stack(files; child_kwargs=(mappedcrs=EPSG(4326),))
stack[:relhum][Lat(Contains(-37), Lon(Contains(144))
```
"""
function stack end

stack(args...; kw...) = GeoStack(args...; kw...)

function filekey(filename::AbstractString)
    ext = splitext(filename)[2]
    if ext in extensions.NCD
        _ncdfilekey(filename)
    else
        throw(ArgumentError("Cannot retreve layer key for $ext files, provide a `keys` argument or filenames as a NamedTuple"))
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
    T = _sourcetype(filename)
    if hasstackfile(T)
        write(filename, _sourcetype(filename), s; kw...)
    else
        # Otherwise write separate arrays
        for key in keys(s)
            fn = joinpath(string(base, "_", key, ext))
            write(fn, _sourcetype(filename), s[key])
        end
    end
end

hasstackfile(T) = false

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
    GeoSeries(filepaths, dims; child=child, kw...)
end

# Support methods

const EXT = Dict(_GRD=>(".grd", ".gri"), _NCD=>(".nc",), _SMAP=>(".h5",))
const REV_EXT = Dict(".grd"=>_GRD, ".gri"=>_GRD, ".nc"=>_NCD, ".h5"=>_SMAP)

_sourcetype(filename::AbstractString) = get(REV_EXT, splitext(filename)[2], _GDAL)

# The the constructor for a geoarray or stack, based on the
# filename extension. GDAL is the fallback for geoarray as it
# handles so many file types.
function _constructor(method::typeof(geoarray), filename; throw=true)
    _, extension = splitext(filename)
    return if extension in EXT[_GRD]
        grdarray
    elseif extension in EXT[_NCD]
        # _check_imported(:NCDatasets, _NCD, extension; throw=throw)
        ncdarray
    elseif extension in EXT[_SMAP]
        # In future we may need to examine the file and check if
        # it's a SMAP file or something else that uses .h5
        throw ? _no_gearray_error(extension) : nothing
    else # GDAL handles too many extensions to list, so just try it and see if it works
        # _check_imported(:ArchGDAL, _GDAL, extension; throw=throw)
        gdalarray
    end
end
_constructor(method::typeof(stack), filename; throw=true) = GeoStack

function _read(f, filename::AbstractString, args...)
    ext = splitext(filename)[2]
    source = get(REV_EXT, ext, _GDAL)
    _read(f, source, filename, args...)
end
_read(f, A::FileArray{X}, args...) where X = _read(f, X, filename(A), args...)

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
