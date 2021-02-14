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
series[Time(Near(DateTime(2001, 1))][:temp][Lat(Between(70, 150)), Lon(Between(-20,20))] |> plot`
```

`GeoSeries` is the only concrete implementation. It includes a `chiltype` field
indicating the constructor used then loading stacks or arrays of any type from disk,
and holds a `kwargs` `NamedTuple` that will be splatted into to the keyword arguments
of the `childtype` constructor. This gives control over the construction of lazy-loaded
files.
"""
abstract type AbstractGeoSeries{T,N,D,A,C} <: AbstractDimensionalArray{T,N,D,A} end

# Interface methods ####################################################

childtype(A::AbstractGeoSeries) = A.childtype
childkwargs(A::AbstractGeoSeries) = A.childkwargs

DD.metadata(A::AbstractGeoSeries) = nothing
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
DD.modify(f, A::AbstractGeoSeries) = rebuild(A, map(child -> modify(f, child), values(A)))

Base.values(A::AbstractGeoSeries) = [A[I] for I in CartesianIndices(A)]

# Array interface methods ##############################################
# Mostly these inherit from AbstractDimensionalArray

@propagate_inbounds function Base.getindex(
    A::AbstractGeoSeries{<:AbstractString}, I::CartesianIndex
)
    childtype(A)(data(A)[I]; refdims=DD.slicedims(A, Tuple(I))[2], A.childkwargs...)
end
@propagate_inbounds function Base.getindex(
    A::AbstractGeoSeries{<:AbstractString}, I::Integer...
)
    childtype(A)(data(A)[I...]; refdims=DD.slicedims(A, I)[2], A.childkwargs...)
end
# Window is passed on to existing MemStacks, as with DiskStacks
@propagate_inbounds function Base.getindex(
    A::AbstractGeoSeries{<:AbstractGeoStack}, I::Integer...
)
    rebuild(data(A)[I...]; A.childkwargs...)
end


"""
    GeoSeries <: AbstractGeoSeries

    GeoSeries(A::AbstractArray{<:AbstractGeoArray}, dims; kw...)
    GeoSeries(A::AbstractArray{<:AbstractGeoStack}, dims; kw...)
    GeoSeries(filenames::AbstractArray{<:AbstractString}, dims; kw...)

Concrete implementation of [`AbstractGeoSeries`](@ref).
Series hold paths to array or stack files, along some dimension(s).

# Keywords

- `refdims`: existing reference 
- `childtype`: type of child objects - an `AbstractGeoSeries` or `AbstractGeoStack`
- `childkwargs`: keyword arguments passed to the child object on construction.
"""
struct GeoSeries{T,N,D,R,A<:AbstractArray{T,N},C,K} <: AbstractGeoSeries{T,N,D,A,C}
    data::A
    dims::D
    refdims::R
    childtype::C
    childkwargs::K
end
function GeoSeries(
    data::Array{T}, dims; refdims=(), childtype=DD.basetypeof(T), childkwargs=()
) where T<:Union{<:AbstractGeoStack,<:AbstractGeoArray}
    GeoSeries(data, DD.formatdims(data, dims), refdims, childtype, childkwargs)
end
function GeoSeries(data, dims; refdims=(), childtype, childkwargs=())
    GeoSeries(data, DD.formatdims(data, dims), refdims, childtype, childkwargs)
end

@inline function DD.rebuild(A::GeoSeries, data, dims::Tuple, refdims, args...)
    GeoSeries(data, dims, refdims, childtype(A), childkwargs(A))
end

@propagate_inbounds function Base.setindex!(A::GeoSeries, x, I::StandardIndices...)
    setindex!(data(A), x, I...)
end
