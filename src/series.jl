"""
An AbstractGeoArray that holds or points to a series of stacks.

This is useful abstraction where data is broken into separate
files accross one or more dimensions, and need to be loaded 
separately. 
"""
abstract type AbstractGeoSeries{T,N,D,C} <: AbstractGeoArray{T,N,D} end

(::Type{T})(data, dims; refdims=(), metadata=Dict(), childtype=GeoArray) where T<:AbstractGeoSeries = 
    T(data, dims, refdims, metadata, childtype)

@inline rebuild(s::T, data, newdims, newrefdims) where T <: AbstractGeoSeries =
    T(data, newdims, newrefdims, metadata(s), childtype(s))

childtype(series::AbstractGeoSeries) = series.childtype

Base.getindex(s::AbstractGeoSeries{<:LazySubArray}, i::Vararg{Integer}) = 
    lsa = parent(s)
    childtype(s)(LazySubArray(lsa[i...], indices(lsa)))
end
Base.getindex(s::AbstractGeoSeries{<:AbstractString}, i::Vararg{Integer}) = 
    childtype(s)(parent(s)[i...])
Base.getindex(s::AbstractGeoSeries{<:AbstractGeoStack}, i::Vararg{Integer}) = 
    parent(s)[i...]

"""
Holds stacks along some dimension(s)

Probably not as useful as format-specific AbstractLazySeries that lazily
loads specific files. 
"""
struct GeoSeries{T,N,D,R,A<:AbstractArray{T,N},M,C} <: AbstractGeoSeries{T,N,D,C}
    data::A
    dims::D
    refdims::R
    metadata::M
    childtype::C
end
