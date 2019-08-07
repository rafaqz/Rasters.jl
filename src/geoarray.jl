struct GeoArray{T,N,D,R,A<:AbstractArray{T,N},Mi,U,Me} <: AbstractGeoArray{T,N,D}
    data::A
    dims::D
    refdims::R
    missingval::Mi
    units::U
    metadata::Me
end
GeoArray(a::AbstractArray{T,N}, dims; refdims=(), missingval=missing, units="", 
         metadata=Dict()) where {T,N} = 
    GeoArray(a, formatdims(a, dims), refdims, missingval, units, metadata)

# Getters
missingval(a::GeoArray) = a.missingval
metadata(a::GeoArray) = a.metadata

# Interfaces
Base.parent(a::GeoArray) = a.data

CoordinateReferenceSystemsBase.crs(a::GeoArray) = get(metadata(a), :crs, nothing)

DimensionalData.dims(a::GeoArray) = a.dims
DimensionalData.refdims(a::GeoArray) = a.refdims
DimensionalData.units(a::GeoArray) = a.units
DimensionalData.name(a::GeoArray) = isnothing(metadata(a)) ? "" : get(metadata(a), :name, "")
DimensionalData.shortname(a::GeoArray) = isnothing(metadata(a)) ? "" : get(metadata(a), :shortname, "")
DimensionalData.rebuild(a::GeoArray, data, dims, refdims, missingval=missingval(a)) =
    GeoArray(data, dims, refdims, missingval, a.units, a.metadata)
