"""
An `AbstractDimensionalArray` that holds or points to a series of stacks.

`AbstractGeoSeries` are a high-level `DimensionalArray`s that hold stacks or 
arrays or the paths they can be loaded from. `GeoSeries` are indexed with dimensions 
as with a `AbstractGeoArray`. This is useful when you have multiple files containing
rasters or stacks of rasters spread over dimensions like time and elevation.
As much as possible, implementations should facilitate loading entire
directories and detecting the dimensions from metadata.

This allows 
```julia
series[Time(Near(DateTime(2001, 1))][:temp][Lat(Between(70, 150)), Lon(Between(-20,20))] |> plot`
```

`GeoSeries` is the only concrete implementation, as it includes a field indicating its
child constructor used if loading stacks or arrays of any type from disk.
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
GeoSeries(data, dims; refdims=(), childtype=GeoStack, kwargs...) =
    GeoSeries(data, formatdims(data, dims), refdims, childtype, kwargs)

@inline rebuild(A::GeoSeries, data, dims::Tuple, refdims, args...) =
    GeoSeries(data, dims, refdims, childtype(A), kwargs(A))

Base.@propagate_inbounds Base.setindex!(A::GeoSeries, x, I::Vararg{<:Union{AbstractArray,Colon,Real}}) =
    setindex!(data(A), x, I...)
