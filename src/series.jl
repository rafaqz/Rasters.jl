"""
`AbstractGeoSeries` are a high-level `DimensionalArray` that hold stacks, arrays, 
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
kwargs(A::AbstractGeoSeries) = A.kwargs
metadata(A::AbstractGeoSeries) = nothing
name(A::AbstractGeoSeries) = ""
label(A::AbstractGeoSeries) = ""

# Array interface methods ##############################################
# Mostly these inherit from AbstractDimensionalArray

Base.getindex(A::AbstractGeoSeries{<:AbstractString}, I::Vararg{<:Integer}) =
    childtype(A)(data(A)[I...]; refdims=slicedims(A, I)[2], A.kwargs...)
# TODO how should window be passed on to existing stacks?
Base.getindex(A::AbstractGeoSeries{<:AbstractGeoStack}, I::Vararg{<:Integer}) =
    rebuild(data(A)[I...]; refdims=slicedims(A, I)[2], A.kwargs...)


"""
    GeoSeries(data::Array{T}, dims; refdims=(), childtype=DD.basetypeof(T), 
              kwargs=()) where T<:Union{<:AbstractGeoStack,<:AbstractGeoArray}
    GeoSeries(data, dims; refdims=(), childtype, kwargs)

Concrete implementation of [`AbstractGeoSeries`](@ref).
Series hold paths to array or stack files, along some dimension(s).
"""
struct GeoSeries{T,N,D,R,A<:AbstractArray{T,N},C,K} <: AbstractGeoSeries{T,N,D,A,C}
    data::A
    dims::D
    refdims::R
    childtype::C
    kwargs::K
end
GeoSeries(data::Array{T}, dims; refdims=(), childtype=DD.basetypeof(T), kwargs=()
         ) where T<:Union{<:AbstractGeoStack,<:AbstractGeoArray} =
    GeoSeries(data, formatdims(data, dims), refdims, childtype, kwargs)
GeoSeries(data, dims; refdims=(), childtype, kwargs=()) =
    GeoSeries(data, formatdims(data, dims), refdims, childtype, kwargs)

@inline rebuild(A::GeoSeries, data, dims::Tuple, refdims, args...) =
    GeoSeries(data, dims, refdims, childtype(A), kwargs(A))

Base.@propagate_inbounds Base.setindex!(A::GeoSeries, x, I::Union{AbstractArray,Colon,Integer}...) =
    setindex!(data(A), x, I...)
