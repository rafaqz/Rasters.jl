"""
Spatial array types that can be indexed using dimensions.
"""
abstract type AbstractGeoArray{T,N,D,A} <: AbstractDimensionalArray{T,N,D,A} end

# Marker singlton for lazy loaded arrays, only used for broadcasting
# Can be removed when DiskArrays.jl is used everywhere
struct LazyArray{T,N} <: AbstractArray{T,N} end


# Interface methods ###########################################################

data(a::AbstractGeoArray) = a.data
dims(a::AbstractGeoArray) = a.dims
refdims(a::AbstractGeoArray) = a.refdims
metadata(a::AbstractGeoArray) = a.metadata
missingval(a::AbstractGeoArray) = a.missingval
window(a::AbstractGeoArray) = a.window
name(a::AbstractGeoArray) = a.name
units(a::AbstractGeoArray) = getmeta(a, :units, "")
label(a::AbstractGeoArray) = string(name(a), " ", units(a))

crs(a::AbstractGeoArray, dim) = crs(dims(a, dim))
crs(dim::AbstractDimension) = crs(metadata(dim))

# Rebuild as GeoArray by default
rebuild(a::AbstractGeoArray, data, dims, refdims) =
    GeoArray(data, dims, refdims, metadata(a), missingval(a), name(a))
rebuild(a::AbstractGeoArray; data=data(a), dims=dims(a), refdims=refdims(a),
        metadata=metadata(a), missingval=missingval(a), name=name(a)) = begin
    GeoArray(data, dims, refdims, metadata, missingval, name)
end


"""
Abstract supertype for all memory-backed GeoArrays where the data is an array.
"""
abstract type MemGeoArray{T,N,D,A} <: AbstractGeoArray{T,N,D,A} end


"""
Abstract supertype for all disk-backed GeoArrays. 
For these the data is lazyily loaded from disk.
"""
abstract type DiskGeoArray{T,N,D,A} <: AbstractGeoArray{T,N,D,A} end

filename(A::DiskGeoArray) = A.filename
Base.size(A::DiskGeoArray) = A.size
window(A::DiskGeoArray) = A.window

Base.write(A::T) where T <: DiskGeoArray = write(filename(A), A)
Base.write(filename::AbstractString, A::T) where T <: DiskGeoArray =
    write(filename, basetypeof(T), A)
Base.write(::Type{T}, A::DiskGeoArray) where T <: DiskGeoArray =
    write(filename(A), T, A)


# Concrete implementation ######################################################

"""
A generic, memory-backed spatial array type.
"""
struct GeoArray{T,N,D<:Tuple,R<:Tuple,A<:AbstractArray{T,N},Me,Mi,Na} <: MemGeoArray{T,N,D,A}
    data::A
    dims::D
    refdims::R
    metadata::Me
    missingval::Mi
    name::Na
end

@inline GeoArray(A::AbstractArray{T,N}, dims; refdims=(), metadata=NamedTuple(),
                 missingval=missing, name=Symbol("")) where {T,N} =
    GeoArray(A, formatdims(A, dims), refdims, metadata, missingval, name)

@inline GeoArray(A::MemGeoArray; data=data(A), dims=dims(A), refdims=refdims(A),
                 metadata=metadata(A), missingval=missingval(A), name=name(A)) =
    GeoArray(data, dims, refdims, metadata, missingval, name)
@inline GeoArray(A::DiskGeoArray; data=data(A), dims=dims(A), refdims=refdims(A),
                 metadata=metadata(A), missingval=missingval(A), name=name(A)) = begin
    _window = maybewindow2indices(A, dims, window(A))
    _dims, _refdims = slicedims(dims, refdims, _window)
    GeoArray(data, _dims, _refdims, metadata, missingval, name)
end

dims(a::GeoArray) = a.dims

Base.@propagate_inbounds Base.setindex!(a::GeoArray, x, I::DimensionalData.StandardIndices) =
    setindex!(data(a), x, I...)

Base.convert(::Type{GeoArray}, array::AbstractGeoArray) = GeoArray(array)

# Manually add broadcast style to GeoArray until all sources are real arrays
# and we have access to type parameter A for all of them.
Base.BroadcastStyle(::Type{<:GeoArray{T,N,D,R,A}}) where {T,N,D,R,A} = begin
    inner_style = typeof(Base.BroadcastStyle(A))
    return DimensionalData.DimensionalStyle{inner_style}()
end

# Utils ########################################################################

@inline getmeta(a::AbstractGeoArray, key, fallback) = getmeta(metadata(a), key, fallback)
@inline getmeta(m::Nothing, key, fallback) = fallback
@inline getmeta(m::Union{NamedTuple,Dict}, key, fallback) = key in keys(m) ?  m[key] : fallback
@inline getmeta(m::Metadata, key, fallback) = getmeta(val(m), key, fallback)
