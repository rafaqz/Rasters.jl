"""
An AbstractGeoArray that holds or points to a series of stacks.

This is useful abstraction where data is broken into separate
files accross one or more dimensions, and need to be loaded 
separately. 
"""
abstract type AbstractGeoSeries{T,N,D,C} <: AbstractGeoArray{T,N,D} end

(::Type{T})(data, dims; refdims=(), metadata=Dict(), childtype=GeoArray, window=()
           ) where T<:AbstractGeoSeries = 
    T(data, dims, refdims, metadata, childtype, window)


childtype(series::AbstractGeoSeries) = series.childtype
window(series::AbstractGeoSeries) = series.window

Base.getindex(s::AbstractGeoSeries{<:AbstractString}, I::Vararg{<:Integer}) = 
    childtype(s)(parent(s)[I...]; refdims=slicedims(s, I)[2], window=window(s))
Base.getindex(series::AbstractGeoSeries{<:AbstractGeoStack}, I::Vararg{<:Integer}) = begin
    stack = parent(series)[I...]
    rebuild(stack; window=window(series))
end

"""
Holds stacks along some dimension(s)

Probably not as useful as format-specific AbstractLazySeries that lazily
loads specific files. 
"""
struct GeoSeries{T,N,D,R,A<:AbstractArray{T,N},M,C,V} <: AbstractGeoSeries{T,N,D,C}
    data::A
    dims::D
    refdims::R
    metadata::M
    childtype::C
    window::V
end

@inline rebuild(s::GeoSeries; parent=parent(s), dims=dims(s), 
                refdims=refdims(s), metadata=metadata(s), window=window(s)) =
    GeoSeries(data, dims, refdims, metadata, childtype, window)
@inline rebuild(s::GeoSeries, data, newdims, newrefdims) =
    GeoSeries(data, newdims, newrefdims, metadata(s), childtype(s), window(s))
