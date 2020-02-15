"""
An AbstractGeoArray that holds or points to a series of stacks.

This is useful abstraction where data is broken into separate
files accross one or more dimensions, and need to be loaded
separately.
"""
abstract type AbstractGeoSeries{T,N,D,A,C} <: AbstractDimensionalArray{T,N,D,A} end

# Interface methods ####################################################

data(A::AbstractGeoSeries) = A.data
childtype(A::AbstractGeoSeries) = A.childtype
window(A::AbstractGeoSeries) = A.window
name(A::AbstractGeoSeries) = ""
label(A::AbstractGeoSeries) = ""

# Array interface methods ##############################################
# Mostly these inherit from AbstractGeoArray

Base.getindex(A::AbstractGeoSeries{<:AbstractString}, I::Vararg{<:Integer}) =
    childtype(A)(data(A)[I...]; refdims=slicedims(A, I)[2], window=window(A))
# TODO how shold window be passed on to existing stacks?
Base.getindex(A::AbstractGeoSeries{<:AbstractGeoStack}, I::Vararg{<:Integer}) =
    rebuild(data(A)[I...]; refdims=slicedims(A, I)[2])


# Concrete implementation ##############################################

"""
Series hold paths to array or stack files, along some dimension(s).
"""
struct GeoSeries{T,N,D,R,A<:AbstractArray{T,N},M,C,V} <: AbstractGeoSeries{T,N,D,A,C}
    data::A
    dims::D
    refdims::R
    metadata::M
    childtype::C
    window::V
end

GeoSeries(data, dims; refdims=(), metadata=Dict(), childtype=GeoStack, window=()) =
    GeoSeries(data, formatdims(data, dims), refdims, metadata, childtype, window)

@inline rebuild(A::GeoSeries, data, newdims, newrefdims, name) =
    GeoSeries(data, newdims, newrefdims, metadata(A), childtype(A), window(A))

Base.@propagate_inbounds Base.setindex!(A::GeoSeries, x, I::Vararg{<:Union{AbstractArray,Colon,Real}}) =
    setindex!(data(A), x, I...)
