"""
    AbstractRasterSeries <: DimensionalData.AbstractDimensionalArray

Abstract supertype for high-level `DimensionalArray` that hold `RasterStacks`, `Rasters`,
or the paths they can be loaded from. `RasterSeries` are indexed with dimensions
as with a `AbstractRaster`. This is useful when you have multiple files containing
rasters or stacks of rasters spread over dimensions like time and elevation.

As much as possible, implementations should facilitate loading entire
directories and detecting the dimensions from metadata.

This allows syntax like for a series of stacks of arrays:

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

"""
    modify(f, series::AbstractRasterSeries)

Apply function `f` to the data of the child object.
If the child is an `AbstractRasterStack` the function will
be passed on to its child `AbstractRaster`s.

`f` must return an idenically sized array.

This method triggers a complete rebuild of all objects,
and disk based objects will be transferred to memory.

This is useful for swapping out array backend for an
entire series to `CuArray` from CUDA.jl to copy data to a GPU,
and potentially other types like `DAarray` from Distributed.jl.
"""
DD.modify(f, A::AbstractRasterSeries) = map(child -> modify(f, child), values(A))

"""
    RasterSeries <: AbstractRasterSeries

    RasterSeries(arrays::AbstractArray{<:AbstractRaster}, dims; kw...)
    RasterSeries(stacks::AbstractArray{<:AbstractRasterStack}, dims; kw...)
    RasterSeries(filepaths::AbstractArray{<:AbstractString}, dims; child, duplicate_first, kw...)
    RasterSeries(dirpath:::AbstractString, dims; ext, child, duplicate_first, kw...)

Concrete implementation of [`AbstractRasterSeries`](@ref).
`RasterSeries` are an array of `Raster` or `RasterStack`, along some dimension(s).

# Arguments

- `dims`: series dimension/s.

# Keywords

- `refdims`: existing reference dimension/s.
- `child`: constructor of child objects for use with filenames are passed in,
    can be `Raster` or `RasterStack`. Defaults to `Raster`.
- `duplicate_first::Bool`: wether to duplicate the dimensions and metadata of the
    first file with all other files. This can save load time with a large
    series where dimensions are essentially identical. `true` by default to improve
    load times. If you need exact metadata, set to `false`.
- `ext`: filename extension such as ".tiff" to find when only a directory path is passed in. 
- `kw`: keywords passed to the child constructor [`Raster`](@ref) or [`RasterStack`](@ref)
    if only file names are passed in.
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
    refdims=(), duplicate_first=true, child=nothing, resize=nothing, kw...
)
    childtype = if isnothing(child)
        eltype(filenames) <: NamedTuple ? RasterStack : Raster
    else
        child
    end
    data = if duplicate_first
        # We assume all dims, metadata and missingvals are the same over the series
        # We just load the first object, and swap in the filenames of the others.
        data1 = if childtype <: AbstractRaster
            childtype(first(filenames); kw...)
        else
            childtype(first(filenames); resize, kw...)
        end
        map(filenames) do fn
            swap_filename(data1, fn)
        end
    else
        # Load everything separately
        if childtype <: AbstractRaster
            [childtype(fn; kw...) for fn in filenames]
        else
            [childtype(fn; resize, kw...) for fn in filenames]
        end
    end
    return RasterSeries(data, DD.format(dims, data); refdims)
end
function RasterSeries(dirpath::AbstractString, dims; ext=nothing, kw...)
    filepaths = filter_ext(dirpath, ext)
    RasterSeries(filepaths, dims; kw...)
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

@deprecate series(args...; kw...) RasterSeries(args...; kw...)

# Swap in the filename/s of an object for another filename, wherever it is.
# This is used to use already loaded metadata of one file with another
# file that is similar or identical besides tha actual raster data.
swap_filename(x, filename) = rebuild(x, data=swap_filename(parent(x), filename))
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
