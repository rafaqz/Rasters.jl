"""
Spacial array types that can be indexed using dimensions.
"""
abstract type AbstractGeoArray{T,N,D} <: AbstractDimensionalArray{T,N,D} end

abstract type AbstractDiskGeoArray{T,N,D} <: AbstractGeoArray{T,N,D} end

@inline refdims(a::AbstractGeoArray) = a.refdims
@inline metadata(a::AbstractGeoArray) = a.metadata
@inline missingval(a::AbstractGeoArray) = a.missingval
@inline rebuild(a::T, data, dims, refdims, missingval=missingval(a)) where T<:AbstractGeoArray = 
    GeoArray(data, dims, refdims, metadata(a), missingval)

units(a::AbstractGeoArray) = getmeta(a, :units, "")  
longname(a::AbstractGeoArray) = getmeta(a, :longname, "")
shortname(a::AbstractGeoArray) = getmeta(a, :shortname, "")
mask(a::AbstractGeoArray) = parent(a) .!= missingval(a)
mask(a::AbstractGeoArray{<:Union{Missing}}) = (!).(ismissing.(parent(a)))

replace_missing(a::AbstractGeoArray, mv) = 
    rebuild(a, replace(a, missingval(a) => mv), dims(a), refdims(a), mv)

CoordinateReferenceSystemsBase.crs(a::AbstractGeoArray) = getmeta(a, :crs, nothing)

Base.parent(a::AbstractGeoArray) = a.data

@mix struct GeoArrayMixin{T,N,D,R,A<:AbstractArray{T,N},Me,Mi}
    data::A
    dims::D
    refdims::R
    metadata::Me
    missingval::Mi
end

"""
A generic, memory-backed spatial array type.
"""
@GeoArrayMixin struct GeoArray{} <: AbstractGeoArray{T,N,D} end

@inline GeoArray(a::A, dims::D, refdims::R, metadata::Me, missingval::Mi
        ) where {A<:AbstractArray{T,N},D,R,Me,Mi} where {T,N} = begin
    dims = formatdims(a, dims)
    GeoArray{T,N,typeof(dims),R,A,Me,Mi}(a, dims, refdims, metadata, missingval)
end

@inline GeoArray(a::AbstractArray{T,N}, dims; 
                 refdims=(), metadata=Dict(), missingval=missing) where {T,N} = 
    GeoArray(a, formatdims(a, dims), refdims, metadata, missingval)
# Move disk backed array to memory
@inline GeoArray(a::AbstractGeoArray) = 
    GeoArray(parent(a), dims(a), refdims(a), metadata(a), missingval(a))

Base.convert(::Type{GeoArray}, array::AbstractGeoArray) = GeoArray(array)

@inline getmeta(a::AbstractGeoArray, key, fallback) = getmeta(metadata(a), key, fallback)
@inline getmeta(m::Nothing, key, fallback) = fallback
@inline getmeta(m::Union{NamedTuple,Dict}, key, fallback) = key in keys(m) ?  m[key] : fallback
