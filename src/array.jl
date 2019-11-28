"""
Spatial array types that can be indexed using dimensions.
"""
abstract type AbstractGeoArray{T,N,D} <: AbstractDimensionalArray{T,N,D} end

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
rebuild(a::AbstractGeoArray; data=parent(a), dims=dims(a), refdims=refdims(a), 
        metadata=metadata(a), missingval=missingval(a), name=name(a)) =
    GeoArray(data, dims, refdims, metadata, missingval, name)

abstract type MemGeoArray{T,N,D} <: AbstractGeoArray{T,N,D} end

abstract type DiskGeoArray{T,N,D} <: AbstractGeoArray{T,N,D} end

filename(A::DiskGeoArray) = A.filename
Base.size(A::DiskGeoArray) = A.size
window(a::DiskGeoArray) = a.window

# Base/Other methods ###########################################################

Base.parent(a::AbstractGeoArray) = a.data


# Concrete implementation ######################################################

"""
A generic, memory-backed spatial array type.
"""
struct GeoArray{T,N,A<:AbstractArray{T,N},D<:Tuple,R<:Tuple,Me,Mi,Na} <: AbstractGeoArray{T,N,D}
    data::A
    dims::D
    refdims::R
    metadata::Me
    missingval::Mi
    name::Na
end

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
boolmask(A::AbstractArray) = boolmask(A, missing)
boolmask(A::AbstractGeoArray) = boolmask(A, missingval(A))
boolmask(A::AbstractGeoArray, missingval) = parent(A) .!== missingval

missingmask(A::AbstractArray) = missingmask(A, missing)
missingmask(A::AbstractGeoArray) = missingmask(A, missingval(A))
missingmask(A::AbstractGeoArray, missingval) = 
    (a -> a === missingval ? missing : false).(parent(A))

"""
    replace_missing(a::AbstractGeoArray, newmissing) 

Replace missing values in the array with a new missing value, also
updating the missingval field.
"""
replace_missing(a::AbstractGeoArray, newmissing) = begin
    data = if ismissing(missingval(a))
        collect(Missings.replace(parent(a), newmissing))
    else
        replace(parent(a), missingval(a) => newmissing)
    end
    rebuild(a; data=data, missingval=newmissing)
end


# Utils ########################################################################

@inline getmeta(a::AbstractGeoArray, key, fallback) = getmeta(metadata(a), key, fallback)
@inline getmeta(m::Nothing, key, fallback) = fallback
@inline getmeta(m::Union{NamedTuple,Dict}, key, fallback) = key in keys(m) ?  m[key] : fallback
@inline getmeta(m::Metadata, key, fallback) = getmeta(val(m), key, fallback) 
