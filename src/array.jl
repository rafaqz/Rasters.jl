"""
Spacial array types that can be indexed using dimensions.
"""
abstract type AbstractGeoArray{T,N,D} <: AbstractDimensionalArray{T,N,D} end

metadata(a::AbstractGeoArray) = a.metadata
missingval(a::AbstractGeoArray) = a.missingval

replace_missing(a::AbstractGeoArray, x) = begin
    newa = replace(a, missing=>x)
    rebuild(newa, parent(newa), dims(newa), refdims(newa), missingval(newa))
end


"""
A generic, memory-backed spacial array type.
"""
struct GeoArray{T,N,D,R,A<:AbstractArray{T,N},Me,Mi} <: AbstractGeoArray{T,N,D}
    data::A
    dims::D
    refdims::R
    metadata::Me
    missingval::Mi
end
GeoArray(a::A, dims::D, refdims::R, metadata::Me, missingval::Mi
        ) where {A<:AbstractArray{T,N},D,R,Me,Mi} where {T,N} = begin
    dims = formatdims(a, dims)
    GeoArray{T,N,typeof(dims),R,A,Me,Mi}(a, dims, refdims, metadata, missingval)
end

GeoArray(a::AbstractArray{T,N}, dims; refdims=(), metadata=Dict(), missingval=missing
        ) where {T,N} = 
    GeoArray(a, formatdims(a, dims), refdims, metadata, missingval)

# Interfaces
Base.parent(a::GeoArray) = a.data
Base.convert(::Type{GeoArray}, a::GeoArray) = a

metadata(a::GeoArray) = a.metadata
refdims(a::GeoArray) = a.refdims
rebuild(a::GeoArray, data, dims, refdims, missingval=missingval(a)) =
    GeoArray(data, dims, refdims, missingval, metadata(a))
units(a::GeoArray) = getmeta(a, :units, "")  
name(a::GeoArray) = getmeta(a, :name, "")
shortname(a::GeoArray) = getmeta(a, :shortname, "")

CoordinateReferenceSystemsBase.crs(a::GeoArray) = get(metadata(a), :crs, nothing)

getmeta(a::AbstractGeoArray, key, fallback) = getmeta(metadata(a), key, fallback)
getmeta(m::Nothing, key, fallback) = fallback
getmeta(m::AbstractMetadata, key, fallback) = get(m, key, fallback)
