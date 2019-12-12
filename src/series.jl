"""
An AbstractGeoArray that holds or points to a series of stacks.

This is useful abstraction where data is broken into separate
files accross one or more dimensions, and need to be loaded 
separately. 
"""
abstract type AbstractGeoSeries{T,N,D,CT} <: AbstractGeoArray{T,N,D} end

# Interface methods ####################################################

childtype(series::AbstractGeoSeries) = series.childtype
childdims(s::AbstractGeoSeries) = s.childdims
window(series::AbstractGeoSeries) = series.window

# Array interface methods ##############################################
# Mostly these inherit from AbstractGeoArray

Base.getindex(s::AbstractGeoSeries{<:AbstractString}, I::Vararg{<:Integer}) = 
    childtype(s)(parent(s)[I...]; dims=childdims(s), refdims=slicedims(s, I)[2], 
                 window=window(s), metadata=metadata(s))
Base.getindex(s::AbstractGeoSeries{<:AbstractGeoStack}, I::Vararg{<:Integer}) = 
    rebuild(parent(s)[I...]; window=window(s), refdims=slicedims(s, I)[2])


# Concrete implementation ##############################################

"""
Series hold paths to array or stack files, along some dimension(s).
"""
struct GeoSeries{T,N,D,R,A<:AbstractArray{T,N},M,CT,CD,V} <: AbstractGeoSeries{T,N,D,CT}
    data::A
    dims::D
    refdims::R
    metadata::M
    childtype::CT
    childdims::CD
    window::V
end

GeoSeries(data, dims; refdims=(), metadata=Dict(), childtype=GeoStack, childdims=(), window=()) = 
    GeoSeries(data, formatdims(data, dims), refdims, metadata, childtype, childdims, window)

@inline rebuild(s::GeoSeries, data, newdims, newrefdims) =
    GeoSeries(data, newdims, newrefdims, metadata(s), childtype(s), childdims(s), window(s))

Base.@propagate_inbounds Base.setindex!(a::GeoSeries, x, I::Vararg{<:Union{AbstractArray,Colon,Real}}) =
    setindex!(parent(a), x, I...)
