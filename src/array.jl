"""
Spacial array types that can be indexed using dimensions.
"""
abstract type AbstractGeoArray{T,N,D} <: AbstractDimensionalArray{T,N,D} end

"""
We specifically don't include A in the mixing: it may not actually be an 
AbstractArray
"""
@premix struct GeoArrayMixin{T,N,D<:Tuple,R<:Tuple,Me,Mi,Na}
    data::A
    dims::D
    refdims::R
    metadata::Me
    missingval::Mi
    name::Na
end

dims(a::AbstractGeoArray) = a.dims
refdims(a::AbstractGeoArray) = a.refdims
metadata(a::AbstractGeoArray) = a.metadata
missingval(a::AbstractGeoArray) = a.missingval
window(a::AbstractGeoArray) = a.window
name(a::AbstractGeoArray) = a.name
units(a::AbstractGeoArray) = getmeta(a, :units, "")  
label(a::AbstractGeoArray) = string(name(a), " ", units(a))

rebuild(a::AbstractGeoArray, parent, dims, refdims) =
    GeoArray(parent, dims, refdims, metadata(a), missingval(a), name(a))
rebuild(a::AbstractGeoArray; parent=parent(a), dims=dims(a), refdims=refdims(a), missingval=missingval(a)) =
    GeoArray(parent, dims, refdims, metadata(a), missingval, name(a))

CoordinateReferenceSystemsBase.crs(a::AbstractGeoArray) = getmeta(a, :crs, nothing)

Base.parent(a::AbstractGeoArray) = a.data

"""
A generic, memory-backed spatial array type.
"""
@GeoArrayMixin struct GeoArray{A<:AbstractArray{T,N}} <: AbstractGeoArray{T,N,D} end

@inline GeoArray(a::A, dims::D, refdims::R, metadata::Me, missingval::Mi, name::Na
        ) where {A<:AbstractArray{T,N},D<:Tuple,R<:Tuple,Me,Mi,Na} where {T,N} = begin
    dims = formatdims(a, dims)
    GeoArray{T,N,typeof(dims),R,Me,Mi,Na,A}(a, dims, refdims, metadata, missingval, name)
end

@inline GeoArray(a::AbstractArray{T,N}, dims; 
                 refdims=(), 
                 metadata=NamedTuple(), 
                 missingval=missing, 
                 name=Symbol("")) where {T,N} = 
    GeoArray(a, formatdims(a, dims), refdims, metadata, missingval, name)
@inline GeoArray(a::AbstractGeoArray) = 
    GeoArray(parent(a), dims(a), refdims(a), metadata(a), missingval(a), name(a))

Base.convert(::Type{GeoArray}, array::AbstractGeoArray) = GeoArray(array)

dims(a::GeoArray) = a.dims

# utils

@inline getmeta(a::AbstractGeoArray, key, fallback) = getmeta(metadata(a), key, fallback)
@inline getmeta(m::Nothing, key, fallback) = fallback
@inline getmeta(m::Union{NamedTuple,Dict}, key, fallback) = key in keys(m) ?  m[key] : fallback

mask(a::AbstractGeoArray) = mask(a, missingval(a))
mask(a::AbstractArray) = mask(a, missing)
mask(a::AbstractGeoArray, missingval) = parent(a) .!= missingval
mask(a::AbstractGeoArray, ::Missing) = .!(ismissing.(parent(a)))

replace_missing(a::AbstractGeoArray, mv) = 
    rebuild(a; parent=replace(a, missingval(a) => mv), missingval=mv)
