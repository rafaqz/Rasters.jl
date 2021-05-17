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
abstract type AbstractGeoSeries{T,N,D,A,C} <: AbstractDimensionalArray{T,N,D,A} end

# Interface methods ####################################################

window(A::AbstractGeoSeries) = A.window
child(A::AbstractGeoSeries) = A.child
childkw(A::AbstractGeoSeries) = A.childkw

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
    x = child(A)(data(A)[I]; refdims=DD.slicedims(A, Tuple(I))[2], A.childkw...)
    _maybeview(x, window(A))
end
@propagate_inbounds function Base.getindex(
    A::AbstractGeoSeries{<:AbstractString}, I::Integer...
)
    x = child(A)(data(A)[I...]; refdims=DD.slicedims(A, I)[2], A.childkw...)
    _maybeview(x, window(A))
end
@propagate_inbounds function Base.getindex(
    A::AbstractGeoSeries{<:Union{<:AbstractGeoArray,<:AbstractGeoStack}}, I::Integer...
)
    _maybeview(data(A)[I...], window(A))
end

_maybeview(A, window::Nothing) = A
_maybeview(A, window) = view(A, window...)


"""
    GeoSeries <: AbstractGeoSeries

    GeoSeries(A::AbstractArray{<:AbstractGeoArray}, dims; kw...)
    GeoSeries(A::AbstractArray{<:AbstractGeoStack}, dims; kw...)
    GeoSeries(filenames::AbstractArray{<:AbstractString}, dims; kw...)

Concrete implementation of [`AbstractGeoSeries`](@ref).
Series hold paths to array or stack files, along some dimension(s).

# Keywords

- `refdims`: existing reference dimension/s 
- `child`: constructor of child objects - `geoarray` or `stack`
"""
struct GeoSeries{T,N,D,R,A<:AbstractArray{T,N},W,C,K} <: AbstractGeoSeries{T,N,D,A,C}
    data::A
    dims::D
    refdims::R
    window::W
    child::C
    childkw::K
end
function GeoSeries(
    data::Array{T}, dims; refdims=(), window=nothing, child=nothing, kw...
) where T<:Union{<:AbstractGeoStack,<:AbstractGeoArray}
    childkw = (; kw...)
    GeoSeries(data, DD.formatdims(data, dims), refdims, window, child, childkw)
end
function GeoSeries(data, dims; refdims=(), window=nothing, child, kw...)
    childkw = (; kw...)
    GeoSeries(data, DD.formatdims(data, dims), refdims, window, child, childkw)
end


@inline function DD.rebuild(
    A::GeoSeries, data, dims::Tuple, refdims=(), name=nothing, 
    window=window(A), child=child(A), childkw=childkw(A)
)
    GeoSeries(data, dims, refdims, window, child, childkw)
end
@inline function DD.rebuild(
    A::GeoSeries; 
    data=data(A), dims=dims(A), refdims=refdims(A), name=nothing, 
    window=window(A), child=child(A), childkw=childkw(A)
)
    GeoSeries(data, dims, refdims, window, child, childkw)
end

@propagate_inbounds function Base.setindex!(A::GeoSeries, x, I::StandardIndices...)
    setindex!(data(A), x, I...)
end
