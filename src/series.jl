"""
Holds or points to a series of stacks over one or more dimensions.

This is useful abstraction where data is broken into separate
files accross one or more dimensions, and need to be loaded 
separately. 
"""
abstract type AbstractGeoSeries{T,N,D} <: AbstractDimensionalArray{T,N,D} end

DimensionalData.dims(a::AbstractGeoSeries) = a.dims

"""
Holds stacks along some dimension(s)
"""
struct GeoSeries{T,N,D,R,A<:AbstractArray{T,N},M} <: AbstractGeoSeries{T,N,D}
    stacks::A
    dims::D
    refdims::R
    metadata::M
end

Base.parent(s::GeoSeries) = s.stacks

DimensionalData.refdims(a::GeoSeries) = a.refdims
DimensionalData.metadata(a::GeoSeries) = a.metadata
DimensionalData.rebuild(s::GeoSeries, data, newdims, newrefdims) =
    GeoSeries(data, newdims, newrefdims, metadata(s))
