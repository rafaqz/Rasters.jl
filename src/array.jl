"""
Spacial array types that can be indexed using dimensions.
"""
abstract type AbstractGeoArray{T,N,D} <: AbstractDimensionalArray{T,N,D} end

DimensionalData.metadata(a::AbstractGeoArray) = a.metadata
missingval(a::AbstractGeoArray) = a.missingval

replace_missing(a::AbstractGeoArray, x) = begin
    newa = replace(a, missing=>x)
    rebuild(newa, parent(newa), dims(newa), refdims(newa), missingval(newa))
end


"""
A generic, memory-backed spacial array type.
"""
struct GeoArray{T,N,D,R,A<:AbstractArray{T,N},Mi,Me} <: AbstractGeoArray{T,N,D}
    data::A
    dims::D
    refdims::R
    missingval::Mi
    metadata::Me
end
GeoArray(a::AbstractArray{T,N}, dims; refdims=(), missingval=missing,
         metadata=Dict()) where {T,N} = 
    GeoArray(a, formatdims(a, dims), refdims, missingval, metadata)

# Interfaces
Base.parent(a::GeoArray) = a.data
Base.convert(::Type{GeoArray}, a::GeoArray) = a

DimensionalData.refdims(a::GeoArray) = a.refdims
DimensionalData.rebuild(a::GeoArray, data, dims, refdims, missingval=missingval(a)) =
    GeoArray(data, dims, refdims, missingval, metadata(a))
DimensionalData.units(a::GeoArray) = getmeta(a, :units, "")  
DimensionalData.name(a::GeoArray) = getmeta(a, :name, "")
DimensionalData.shortname(a::GeoArray) = getmeta(a, :shortname, "")

CoordinateReferenceSystemsBase.crs(a::GeoArray) = get(metadata(a), :crs, nothing)

getmeta(a::AbstractGeoArray, key, fallback) = getmeta(metadata(a), key, fallback)
getmeta(m::Nothing, key, fallback) = fallback
getmeta(m::AbstractMetadata, key, fallback) = get(m, key, fallback)
