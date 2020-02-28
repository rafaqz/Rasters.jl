"""
Spatial array types that can be indexed using dimensions.
"""
abstract type AbstractGeoArray{T,N,D,A} <: AbstractDimensionalArray{T,N,D,A} end

# Marker singlton for lazy loaded arrays, only used for broadcasting
# Can be removed when DiskArrays.jl is used everywhere
struct LazyArray{T,N} <: AbstractArray{T,N} end


# Interface methods ###########################################################

data(A::AbstractGeoArray) = A.data
dims(A::AbstractGeoArray) = A.dims
refdims(A::AbstractGeoArray) = A.refdims
metadata(A::AbstractGeoArray) = A.metadata
missingval(A::AbstractGeoArray) = A.missingval
window(A::AbstractGeoArray) = A.window
name(A::AbstractGeoArray) = A.name
units(A::AbstractGeoArray) = getmeta(A, :units, "")
label(A::AbstractGeoArray) = string(name(A), " ", units(A))

# TODO fix this
crs(A::AbstractGeoArray, dim) = crs(dims(A, dim))
crs(dim::Dimension) = get(metadata(dim), :crs, nothing)

# Rebuild all types of AbstractGeoArray as GeoArray
rebuild(A::AbstractGeoArray, data, dims::Tuple, refdims, name=name(A)) =
    GeoArray(data, dims, refdims, name, metadata(A), missingval(A))
rebuild(A::AbstractGeoArray; data=data(A), dims=dims(A), refdims=refdims(A),
        name=name(A), metadata=metadata(A), missingval=missingval(A)) =
    GeoArray(data, dims, refdims, name, metadata, missingval)


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
struct GeoArray{T,N,D<:Tuple,R<:Tuple,A<:AbstractArray{T,N},Na<:AbstractString,Me,Mi} <: MemGeoArray{T,N,D,A}
    data::A
    dims::D
    refdims::R
    name::Na
    metadata::Me
    missingval::Mi
end

@inline GeoArray(A::AbstractArray{T,N}, dims;
                 refdims=(),
                 name="",
                 metadata=NamedTuple(),
                 missingval=missing,
                ) where {T,N} =
    GeoArray(A, formatdims(A, dims), refdims, name, metadata, missingval)

@inline GeoArray(A::MemGeoArray; data=data(A), dims=dims(A), refdims=refdims(A),
                 name=name(A), metadata=metadata(A), missingval=missingval(A)) =
    GeoArray(data, dims, refdims, name, metadata, missingval)
@inline GeoArray(A::DiskGeoArray; data=data(A), dims=dims(A), refdims=refdims(A),
                 name=name(A), metadata=metadata(A), missingval=missingval(A)) = begin
    _window = maybewindow2indices(A, dims, window(A))
    _dims, _refdims = slicedims(dims, refdims, _window)
    GeoArray(data, _dims, _refdims, name, metadata, missingval)
end

dims(A::GeoArray) = A.dims

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
