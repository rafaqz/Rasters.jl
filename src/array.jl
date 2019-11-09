"""
Spatial array types that can be indexed using dimensions.
"""
abstract type AbstractGeoArray{T,N,D} <: AbstractDimensionalArray{T,N,D} end

"""
Common fields for AbstractGeoArray. We explicitly *don't* include the 
type parameter A in the mixin: it may not actually be an AbstractArray. 
Implementations must specify the A parameter and decide its type.
"""
@premix struct GeoArrayMixin{T,N,D<:Tuple,R<:Tuple,Me,Mi,Na}
    data::A
    dims::D
    refdims::R
    metadata::Me
    missingval::Mi
    name::Na
end

# Interface methods ###########################################################

dims(a::AbstractGeoArray) = a.dims
refdims(a::AbstractGeoArray) = a.refdims
metadata(a::AbstractGeoArray) = a.metadata
missingval(a::AbstractGeoArray) = a.missingval
window(a::AbstractGeoArray) = a.window
name(a::AbstractGeoArray) = a.name
units(a::AbstractGeoArray) = getmeta(a, :units, "")  
label(a::AbstractGeoArray) = string(name(a), " ", units(a))

rebuild(a::AbstractGeoArray, data, dims, refdims) =
    GeoArray(data, dims, refdims, metadata(a), missingval(a), name(a))
rebuild(a::AbstractGeoArray; data=parent(a), dims=dims(a), refdims=refdims(a), missingval=missingval(a)) =
    GeoArray(data, dims, refdims, metadata(a), missingval, name(a))


# Base/Other methods ###########################################################

crs(a::AbstractGeoArray) = getmeta(a, :crs, nothing)

Base.parent(a::AbstractGeoArray) = a.data


# Concrete implementation ######################################################

"""
A generic, memory-backed spatial array type.
"""
@GeoArrayMixin struct GeoArray{A<:AbstractArray{T,N}} <: AbstractGeoArray{T,N,D} end

@inline GeoArray(a::AbstractArray{T,N}, dims; refdims=(), metadata=NamedTuple(), 
                 missingval=missing, name=Symbol("")) where {T,N} = 
    GeoArray(a, formatdims(a, dims), refdims, metadata, missingval, name)

@inline GeoArray(A::AbstractGeoArray; data=parent(A), dims=dims(A), refdims=refdims(A), 
                 metadata=metadata(A), missingval=missingval(A), name=name(A)) =
    GeoArray(data, dims, refdims, metadata, missingval, name)

dims(a::GeoArray) = a.dims

Base.@propagate_inbounds Base.setindex!(a::GeoArray, x, I::Vararg{<:Union{AbstractArray,Colon,Real}}) =
    setindex!(parent(a), x, I...)

Base.convert(::Type{GeoArray}, array::AbstractGeoArray) = GeoArray(array)


# Helper methods ##############################################################

mask(a::AbstractGeoArray) = mask(a, missingval(a))
mask(a::AbstractArray) = mask(a, missing)
mask(a::AbstractGeoArray, missingval) = parent(a) .!== missingval
mask(a::AbstractGeoArray, ::Missing) = .!(ismissing.(parent(a)))

"""
    replace_missing(a::AbstractGeoArray, newmissing) 

Replace missing values in the array with a new missing value, also
updating the missingval field.
"""
replace_missing(a::AbstractGeoArray, newmissing) = 
    rebuild(a; data=replace(a, missingval(a) => newmissing), missingval=newmissing)


# Utils ########################################################################

@inline getmeta(a::AbstractGeoArray, key, fallback) = getmeta(metadata(a), key, fallback)
@inline getmeta(m::Nothing, key, fallback) = fallback
@inline getmeta(m::Union{NamedTuple,Dict}, key, fallback) = key in keys(m) ?  m[key] : fallback
