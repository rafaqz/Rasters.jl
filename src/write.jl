
"""
    Base.write(filename::AbstractString, A::AbstractRaster; [source], kw...)

Write an [`AbstractRaster`](@ref) to file, guessing the backend from the
file extension or using the `source` keyword.

# Keywords

$FORCE_KEYWORD
$SOURCE_KEYWORD

Other keyword arguments are passed to the `write` method for the backend.
"""
function Base.write( filename::AbstractString, A::AbstractRaster; 
    source=_sourcetype(filename), 
    kw...
)
    source = _sourcetype(source)
    write(filename, source, A; kw...)
end
Base.write(A::AbstractRaster; kw...) = write(filename(A), A; kw...)
# Fallback
function Base.write(
    filename::AbstractString, source::Source, A::Union{AbstractRaster,AbstractRasterStack}; kw...
)
    missing_package = SOURCE2PACKAGENAME[source]
    error("Missing package extension for $source. Run `using $missing_package` before useing `write` for this file.")
end

"""
    Base.write(filename::AbstractString, s::AbstractRasterStack; kw...)

Write any [`AbstractRasterStack`](@ref) to one or multiple files,
depending on the backend. Backend is guessed from the filename extension
or forced with the `source` keyword.

## Keywords

$EXT_KEYWORD
$FORCE_KEYWORD
$SOURCE_KEYWORD
$SUFFIX_KEYWORD
$VERBOSE_KEYWORD

Other keyword arguments are passed to the `write` method for the backend.

If the source can't be saved as a stack-like object, individual array layers will be saved.
"""
function Base.write(path::AbstractString, s::AbstractRasterStack;
    suffix=nothing,
    ext=nothing,
    source=_sourcetype(path, ext),
    verbose=true,
    kw...
)
    source = _sourcetype(source)
    if haslayers(source)
        write(path, source, s; kw...)
    else
        # Otherwise write separate files for each layer
        if isnothing(ext)
            base, ext = splitext(path)
        else
            base, _ = splitext(path)
        end
        suffix1 = if isnothing(suffix)
            divider = Sys.iswindows() ? '\\' : '/'
            # Add an underscore to the key if there is a file name already
            spacer = last(path) == divider ? "" : "_"
            map(k -> string(spacer, k), keys(s))
        else
            suffix
        end
        if verbose
            @warn string("Cannot write complete stacks to \"", ext, "\", writing layers as individual files")
        end
        map(keys(s)) do key
            fn = string(base, suffix1, ext)
            write(fn, source, s[key]; kw...)
        end |> NamedTuple{keys(s)}
    end
end

"""
    Base.write(filepath::AbstractString, s::AbstractRasterSeries; kw...)

Write any [`AbstractRasterSeries`](@ref) to multiple files,
guessing the backend from the file extension.

The lookup values of the series will be appended to the filepath (before the extension),
separated by underscores.

## Keywords

$EXT_KEYWORD
$FORCE_KEYWORD
$SOURCE_KEYWORD
$VERBOSE_KEYWORD

See specific docs for `write` for eg. gdal or netcdf keywords, or for `RasterStack`. 
All keywords are passed through to these `Raster` and `RasterStack` methods.
"""
function Base.write(path::AbstractString, A::AbstractRasterSeries;
    ext=nothing,
    source=_sourcetype(path, ext),
    verbose=true,
    kw...
)
    source = _sourcetype(source)
    base, path_ext = splitext(path)
    ext = isnothing(ext) ? path_ext : ext
    map(A, DimPoints(A)) do raster, p
        lookupstring = join(map(string, p), "_")
        written_paths = if raster isa RasterStack && !haslayers(source)
            stack_dir = joinpath(dirname(base), lookupstring)
            mkpath(stack_dir)
            stack_filepath = joinpath(stack_dir, basename(path))
            write(stack_filepath, raster; source, verbose, ext, kw...)
        else
            dirpath = dirname(base)
            filename = basename(base) == "" ? lookupstring : join((basename(base), lookupstring), '_')
            slice_filename = joinpath(dirpath, string(filename, ext))
            write(slice_filename, raster; source, verbose, kw...)
        end
        # Don't print warnings for every slice of the series
        verbose = false
        written_paths
    end
end

# Trait for source data that has stack layers
haslayers(T) = false

#  This is used in source `write` methods
function check_can_write(filename, force)
    if !check_can_write(Bool, filename, force)
        throw(ArgumentError("filename already exists at $filename. use the keyword `force=true` to write anyway"))
    end
    return true
end
check_can_write(::Type{Bool}, filename, force) = (force || !isfile(filename))
