"""
Spatial array types that can be indexed using dimensions.
"""
abstract type AbstractGeoArray{T,N,D} <: AbstractDimensionalArray{T,N,D} end

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

abstract type MemGeoArray{T,N,D} <: AbstractGeoArray{T,N,D} end

abstract type DiskGeoArray{T,N,D} <: AbstractGeoArray{T,N,D} end

filename(A::DiskGeoArray) = A.filename
Base.size(A::DiskGeoArray) = A.size
window(A::DiskGeoArray) = A.window

Base.write(A::T) where T <: DiskGeoArray = write(filename(A), A)
Base.write(filename::AbstractString, A::T) where T <: DiskGeoArray =
    write(filename, basetypeof(T), A)
Base.write(::Type{T}, A::DiskGeoArray) where T <: DiskGeoArray =
    write(filename(A), T, A)

# Base/Other methods ###########################################################


# Concrete implementation ######################################################

"""
A generic, memory-backed spatial array type.
"""
struct GeoArray{T,N,D<:Tuple,R<:Tuple,A<:AbstractArray{T,N},Me,Mi,Na} <: MemGeoArray{T,N,D}
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


# Helper methods ##############################################################
boolmask(A::AbstractArray) = boolmask(A, missing)
boolmask(A::AbstractGeoArray) =
    rebuild(A; data=boolmask(A, missingval(A)), missingval=false, name="Boolean mask")
boolmask(A::AbstractGeoArray, missingval::Missing) =
    (x -> !ismissing(x)).(data(A))
boolmask(A::AbstractGeoArray, missingval) =
    (x -> !isapprox(x, missingval)).(data(A))

missingmask(A::AbstractArray) = missingmask(A, missing)
missingmask(A::AbstractGeoArray) =
    rebuild(A; data=missingmask(A, missingval(A)), missingval=missing, name="Missing mask")
missingmask(A::AbstractGeoArray, missingval::Missing) =
    (a -> ismissing(a) ? missing : true).(data(A))
missingmask(A::AbstractGeoArray, missingval) =
    (a -> isapprox(a, missingval) ? missing : true).(data(A))


"""
    replace_missing(a::AbstractGeoArray, newmissing)

Replace missing values in the array with a new missing value, also
updating the missingval field.
"""
replace_missing(a::AbstractGeoArray, newmissing=missing) = begin
    newdata = if ismissing(missingval(a))
        collect(Missings.replace(data(a), newmissing))
    else
        replace(data(a), missingval(a) => newmissing)
    end
    rebuild(a; data=newdata, missingval=newmissing)
end

# Utils ########################################################################

@inline getmeta(a::AbstractGeoArray, key, fallback) = getmeta(metadata(a), key, fallback)
@inline getmeta(m::Nothing, key, fallback) = fallback
@inline getmeta(m::Union{NamedTuple,Dict}, key, fallback) = key in keys(m) ?  m[key] : fallback
@inline getmeta(m::Metadata, key, fallback) = getmeta(val(m), key, fallback)
