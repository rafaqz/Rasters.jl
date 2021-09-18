"""
    AbstractGeoSeries <: DimensionalData.AbstractDimensionalArray

Abstract supertype for high-level `DimensionalArray` that hold stacks, arrays,
or the paths they can be loaded from. `GeoSeries` are indexed with dimensions
as with a `AbstractGeoArray`. This is useful when you have multiple files containing
rasters or stacks of rasters spread over dimensions like time and elevation.
As much as possible, implementations should facilitate loading entire
directories and detecting the dimensions from metadata.

This allows syntax like:

```julia
series[Time(Near(DateTime(2001, 1))][:temp][Y(Between(70, 150)), X(Between(-20,20))] |> plot`
```

`GeoSeries` is the only concrete implementation. It includes a `chiltype` field
indicating the constructor used then loading stacks or arrays of any type from disk,
and holds a `kwargs` `NamedTuple` that will be splatted into to the keyword arguments
of the `child` constructor. This gives control over the construction of lazy-loaded
files.
"""
abstract type AbstractGeoSeries{T,N,D,A} <: AbstractDimensionalArray{T,N,D,A} end

# Interface methods ####################################################

DD.metadata(A::AbstractGeoSeries) = NoMetadata()
DD.name(A::AbstractGeoSeries) = NoName()
DD.label(A::AbstractGeoSeries) = ""

"""
    modify(f, series::AbstractGeoSeries)

Apply function `f` to the data of the child object.
If the child is an `AbstractGeoStack` the function will
be passed on to its child `AbstractGeoArray`s.

`f` must return an idenically sized array.

This method triggers a complete rebuild of all objects,
and disk based objects will be transferred to memory.

This is useful for swapping out array backend for an
entire series to `CuArray` from CUDA.jl to copy data to a GPU,
and potentially other types like `DAarray` from Distributed.jl.
"""
DD.modify(f, A::AbstractGeoSeries) = map(child -> modify(f, child), values(A))

"""
    GeoSeries <: AbstractGeoSeries

    GeoSeries(A::AbstractArray{<:AbstractGeoArray}, dims; kw...)
    GeoSeries(A::AbstractArray{<:AbstractGeoStack}, dims; kw...)
    GeoSeries(filenames::AbstractArray{<:AbstractString}, dims; kw...)

Concrete implementation of [`AbstractGeoSeries`](@ref).
Series hold GeoArray or GeoStack files, along some dimension(s).

# Keywords

- `dims` known dimensions. These are usually read from the first file in the series
    and are assumed to be _the same for all stacks/arrays in the series_.
- `refdims`: existing reference dimension/s 
- `child`: constructor of child objects - `GeoArray` or `stack`
"""
struct GeoSeries{T,N,D,R,A<:AbstractArray{T,N}} <: AbstractGeoSeries{T,N,D,A}
    data::A
    dims::D
    refdims::R
end
function GeoSeries(data::AbstractArray{<:Union{AbstractGeoStack,AbstractGeoArray}}, dims; 
    refdims=()
)
    GeoSeries(data, DD.formatdims(data, dims), refdims)
end
function GeoSeries(filenames::AbstractArray{<:Union{AbstractString,NamedTuple}}, dims; 
    duplicate_first=true, child=nothing, resize=nothing, kw...
)
    childtype = if isnothing(child)
        eltype(filenames) <: NamedTuple ? GeoStack : GeoArray
    else
        child
    end
    data = if duplicate_first
        # We assume all dims, metadata and missingvals are the same over the series
        # We just load the first object, and swap in the filenames of the others.
        data1 = if childtype <: AbstractGeoArray
            childtype(first(filenames); kw...)
        else
            childtype(first(filenames); resize, kw...)
        end
        map(filenames) do fn
            swap_filename(data1, fn)
        end
    else
        # Load everything separately
        if childtype <: AbstractGeoArray
            [childtype(fn; kw...) for fn in filenames]
        else
            [childtype(fn; resize, kw...) for fn in filenames]
        end
    end
    return GeoSeries(data, DD.formatdims(data, dims))
end

function GeoSeries(dirpath::AbstractString, dims=(Dim{:series}(),); ext=nothing, child=GeoArray, kw...)
    filepaths = filter_ext(dirpath, ext)
    GeoSeries(filepaths, dims; child=child, kw...)
end

@inline function DD.rebuild(
    A::GeoSeries, data, dims::Tuple, refdims=(), name=nothing, metadata=nothing,
)
    GeoSeries(data, dims, refdims)
end
@inline function DD.rebuild(
    A::GeoSeries; 
    data=data(A), dims=dims(A), refdims=refdims(A), name=nothing, metadata=nothing,
)
    GeoSeries(data, dims, refdims)
end

@deprecate series(args...; kw...) GeoSeries(args...; kw...)

swap_filename(x, filename) = rebuild(x, data=swap_filename(data(x), filename))
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
