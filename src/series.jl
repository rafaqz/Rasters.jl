"""
    AbstractRasterSeries <: DimensionalData.AbstractDimensionalArray

Abstract supertype for high-level `DimensionalArray` that hold `RasterStacks`, `Rasters`,
or the paths they can be loaded from. `RasterSeries` are indexed with dimensions
as with a `AbstractRaster`. This is useful when you have multiple files containing
rasters or stacks of rasters spread over dimensions like time and elevation.

As much as possible, implementations should facilitate loading entire
directories and detecting the dimensions from metadata.

This allows syntax like below for a series of stacks of arrays:

```julia
RasterSeries[Time(Near(DateTime(2001, 1))][:temp][Y(Between(70, 150)), X(Between(-20,20))] |> plot`
```

[`RasterSeries`](@ref) is the concrete implementation.
"""
abstract type AbstractRasterSeries{T,N,D,A} <: AbstractDimArray{T,N,D,A} end

# Interface methods ####################################################

DD.metadata(A::AbstractRasterSeries) = NoMetadata()
DD.name(A::AbstractRasterSeries) = NoName()
DD.label(A::AbstractRasterSeries) = ""
isdisk(A::AbstractRasterSeries) = any(isdisk, A)

"""
    modify(f, series::AbstractRasterSeries)

Apply function `f` to the data of the child object.
If the child is an `AbstractRasterStack` the function will
be passed on to its child `AbstractRaster`s.

`f` must return an identically sized array.

This method triggers a complete rebuild of all objects,
and disk based objects will be transferred to memory.

An example of the usefulnesss of this is for swapping out array backend
for an entire series to `CuArray` from CUDA.jl to copy data to a GPU.
"""
DD.modify(f, A::AbstractRasterSeries) = map(child -> modify(f, child), values(A))

"""
    RasterSeries <: AbstractRasterSeries

    RasterSeries(rasters::AbstractArray{<:AbstractRaster}, dims; [refdims])
    RasterSeries(stacks::AbstractArray{<:AbstractRasterStack}, dims; [refdims]) 

    RasterSeries(paths::AbstractArray{<:AbstractString}, dims; child, duplicate_first, kw...)
    RasterSeries(path:::AbstractString, dims; ext, separator, child, duplicate_first, kw...)

Concrete implementation of [`AbstractRasterSeries`](@ref).

A `RasterSeries` is an array of `Raster`s or `RasterStack`s, along some dimension(s).

Existing `Raster` `RasterStack` can be wrapped in a `RasterSeries`, or new files 
can be loaded from an array of `String` or from a single `String`.

A single `String` can refer to a whole directory, or the name of a series of files in a directory,
sharing a common stem. The differnce between the filenames can be used as the lookup for the
series. 

For example, with some tifs at these paths : 

```
"series_dir/myseries_2001-01-01T00:00:00.tif"
"series_dir/myseries_2002-01-01T00:00:00.tif"
```

We can load a `RasterSeries` with a `DateTime` lookup:

```julia
julia> ser = RasterSeries("series_dir/myseries.tif", Ti(DateTime))
2-element RasterSeries{Raster,1} with dimensions: 
  Ti Sampled{DateTime} DateTime[DateTime("2001-01-01T00:00:00"), DateTime("2002-01-01T00:00:00")] ForwardOrdered Irregular Points
```

The `DateTime` suffix is parsed from the filenames. Using `Ti(Int)` would try to parse integers instead.

Just using the directory will also work, unless there are other files mixed in it:

```julia
julia> ser = RasterSeries("series_dir", Ti(DateTime))
2-element RasterSeries{Raster,1} with dimensions: 
  Ti Sampled{DateTime} DateTime[DateTime("2001-01-01T00:00:00"), DateTime("2002-01-01T00:00:00")] ForwardOrdered Irregular Points
```

# Arguments

- `dims`: series dimension/s.

# Keywords

When loading a series from a Vector of `String` paths or a single `String` path:
- `child`: constructor of child objects for use when filenames are passed in,
    can be `Raster` or `RasterStack`. Defaults to `Raster`.
- `duplicate_first::Bool`: wether to duplicate the dimensions and metadata of the
    first file with all other files. This can save load time with a large
    series where dimensions are identical. `false` by default.
$LAZY_KEYWORD
- `kw`: keywords passed to the child constructor [`Raster`](@ref) or [`RasterStack`](@ref).

When loading a series from a single `String` path:

- `separator`: separator used to split lookup elements from the rest of a filename. '_' by default.


Others:
- `refdims`: existing reference dimension/s, normally not required.
"""
struct RasterSeries{T,N,D,R,A<:AbstractArray{T,N}} <: AbstractRasterSeries{T,N,D,A}
    data::A
    dims::D
    refdims::R
end
function RasterSeries(data::AbstractArray{<:Union{AbstractRasterStack,AbstractRaster}}, dims;
    refdims=()
)
    RasterSeries(data, DD.format(dims, data), refdims)
end
function RasterSeries(filenames::NamedTuple{K}, dims; kw...) where K
    RasterSeries(map((fns...) -> NamedTuple{K}(fns), values(filenames)...), dims; kw...)
end
function RasterSeries(filenames::AbstractArray{<:Union{AbstractString,NamedTuple}}, dims;
    refdims=(), 
    lazy=false, 
    duplicate_first=false, 
    child=nokw, 
    resize=nokw, 
    kw...
)
    childtype = if child isa NoKW
        eltype(filenames) <: NamedTuple ? RasterStack : Raster
    else
        child
    end
    data = if lazy && duplicate_first
        # We assume all dims, metadata and missingvals are the same over the series
        # We just load the first object, and swap in the filenames of the others.
        data1 = if childtype <: AbstractRaster
            childtype(first(filenames); lazy, kw...)
        else
            childtype(first(filenames); lazy, resize, kw...)
        end
        map(filenames) do fn
            swap_filename(data1, fn)
        end
    else
        # Load everything separately
        if childtype <: AbstractRaster
            [childtype(fn; lazy, kw...) for fn in filenames]
        else
            [childtype(fn; resize, lazy, kw...) for fn in filenames]
        end
    end
    return RasterSeries(data, DD.format(dims, data); refdims)
end
function RasterSeries(path::AbstractString, dims; 
    refdims=(), 
    ext=nothing, 
    separator='_', 
    kw...
)
    if isdir(path)
        dirpath = path
        filepaths = filter_ext(dirpath, ext)
        length(filepaths) > 0 || error("No $(isnothing(ext) ? "" : ext) files found in \"$path\" dir")
        full_filename, _ = splitext(basename(first(filepaths)))
        common_filename = join(split(full_filename, separator)[1:end-1])
    else
        dirpath = dirname(path)
        common_filename, path_ext = splitext(basename(path))
        ext = (isnothing(ext) && path_ext != "") ? path_ext : ext
        filepaths = filter(filter_ext(dirpath, ext)) do fp
            basename(fp)[1:length(common_filename)] == common_filename
        end
        length(filepaths) > 0 || error("No $(isnothing(ext) ? "" : ext) files found matching \"$common_filename\" stem and \"$ext\" extension")
    end
    basenames = map(fp -> basename(splitext(fp)[1]), filepaths)
    v = val(dims)
    # Try to get values of the wrapped type from the filenames
    if dims isa Dimension && val(dims) isa Union{Type,AutoVal{<:Type}}
        T = val(dims) isa Type ? v : val(v)
        index_strings = map(basenames) do n
            strip(n[length(common_filename)+1:end], separator)
        end
        index = map(index_strings) do s
            try
                parse(T, s) 
            catch
                error("Could not parse filename segment $s as $T")
            end
        end
        dims = if v isa Type
            (basetypeof(dims)(index),)
        else
            (basetypeof(dims)(index; v.kw...),)
        end
    end
    RasterSeries(filepaths, DD.format(dims, filepaths); refdims, kw...)
end

@inline function DD.rebuild(
    A::RasterSeries, data, dims::Tuple, refdims=(), name=nothing, metadata=nothing,
)
    RasterSeries(data, dims, refdims)
end
@inline function DD.rebuild(
    A::RasterSeries;
    data=parent(A), dims=dims(A), refdims=refdims(A), name=nothing, metadata=nothing,
)
    RasterSeries(data, dims, refdims)
end

function Base.map(f, series::RasterSeries)
    vals = map(f, parent(series))
    if eltype(vals) <: Union{AbstractRaster,AbstractRasterStack}
        return rebuild(series, vals)
    else
        return Raster(vals, dims(series))
    end
end

# Swap in the filename/s of an object for another filename, wherever it is.
# This is used to use already loaded metadata of one file with another
# file that is similar or identical besides tha actual raster data.
swap_filename(x, filename) = rebuild(x, data=swap_filename(parent(x), filename))
swap_filename(x::Raster, filename::AbstractString) = rebuild(x, data=swap_filename(parent(x), filename))
swap_filename(x::NamedTuple, filenames::NamedTuple) = map(swap_filename, x, filenames)
swap_filename(x::FileStack, filename::AbstractString) = @set x.filename = filename
swap_filename(x::FileArray, filename::AbstractString) = @set x.filename = filename
function swap_filename(x::AbstractArray, filename::AbstractString)
    # The `FileArray` is wrapped, so use Flatten.jl to update it wherever it is
    ignore = Union{Dict,Set,Base.MultiplicativeInverses.SignedMultiplicativeInverse}
    Flatten.modify(x, FileArray, ignore) do fa

        @set fa.filename = filename
    end
end
