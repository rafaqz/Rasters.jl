"""
`AbstractGeoArray` wraps an array (or location of an array) and metadata
about its contents. It may be memory (`GeoArray`) or disk-backed (`NCarray`,
`GDAlarray`).

`AbstractGeoArray`s inherit from `AbstractDimensionalArray` from
DimensionalData. They can be indexed as regular Julia arrays or with
DimensionalData.jl dimensions. They will plot as a heatmap in Plots.jl with correct
coordinates and labels, even after slicing with `getindex` or `view`. `getindex`
on a `AbstractGeoArray` will always return a standard `GeoArray`.

 In addition to DimensionalArray behaviour, these have
`metadata` and `missingval` fields
"""
abstract type AbstractGeoArray{T,N,D,A} <: AbstractDimensionalArray{T,N,D,A} end

# Marker singleton for lazy loaded arrays, only used for broadcasting.
# Can be removed when DiskArrays.jl is used everywhere
struct LazyArray{T,N} <: AbstractArray{T,N} end


# Interface methods ###########################################################

data(A::AbstractGeoArray) = A.data
dims(A::AbstractGeoArray) = A.dims
refdims(A::AbstractGeoArray) = A.refdims
name(A::AbstractGeoArray) = A.name
metadata(A::AbstractGeoArray) = A.metadata
missingval(A::AbstractGeoArray) = A.missingval
window(A::AbstractGeoArray) = A.window
units(A::AbstractGeoArray) = getmeta(A, :units, "")
label(A::AbstractGeoArray) = string(name(A), " ", units(A))

"""
    crs(A::AbstractGeoArray)
Get the coordinate reference system of the array.
"""
crs(A::AbstractGeoArray) =
    if hasdim(A, Lat)
        crs(dims(A, Lat))
    elseif hasdim(A, Lon)
        crs(dims(A, Lon))
    else
        error("No Lat or Lon dimension, crs not available")
    end
"""
    crs(A::AbstractGeoArray)
Get the coordinate reference system of a Dimension.
"""
crs(dim::Dimension) = crs(mode(dim), dim)
"""
    usercrs(A::AbstractGeoArray)
Get the coordinate reference system of the array.
"""
usercrs(A::AbstractGeoArray) =
    if hasdim(A, Lat)
        usercrs(dims(A, Lat))
    elseif hasdim(A, Lon)
        usercrs(dims(A, Lon))
    else
        error("No Lat or Lon dimension, usercrs not available")
    end
"""
    usercrs(A::AbstractGeoArray)
Get the user input coordinate reference system.
"""
usercrs(dim::Dimension) = usercrs(mode(dim), dim)

# Rebuild all types of AbstractGeoArray as GeoArray
rebuild(A::AbstractGeoArray, data, dims::Tuple, refdims, name, metadata, missingval=missingval(A)) =
    GeoArray(data, dims, refdims, name, metadata, missingval)
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
A generic, memory-backed spatial array type. All [`AbstractGeoArray`](@ref) are
converted to GeoArray when indexed or otherwise transformed.
"""
struct GeoArray{T,N,D<:Tuple,R<:Tuple,A<:AbstractArray{T,N},Na<:AbstractString,Me,Mi} <: MemGeoArray{T,N,D,A}
    data::A
    dims::D
    refdims::R
    name::Na
    metadata::Me
    missingval::Mi
end
"""
    GeoArray(A::AbstractArray{T,N}, dims::Tuple;
             refdims=(), name="", metadata=nothing, missingval=missing)

Construct a [`GeoArray`](@ref) from an `AbstractArray`, a `Tuple` of
`Dimension` and keyword arguments.
"""
@inline GeoArray(A::AbstractArray{T,N}, dims::Tuple;
                 refdims=(), name="", metadata=nothing, missingval=missing,
                ) where {T,N} =
    GeoArray(A, formatdims(A, dims), refdims, name, metadata, missingval)

"""
    GeoArray(A::AbstractGeoArray; [data=data(A), dims=dims(A), refdims=refdims(A),
             name=name(A), metadata=metadata(A), missingval=missingval(A)]) =

Construct a [`GeoArray`](@ref) from another [`AbstractGeoArray`](@ref), and
keyword arguments.
"""
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

Base.@propagate_inbounds Base.setindex!(A::GeoArray, x, I::Vararg{DimensionalData.StandardIndices}) =
    setindex!(data(A), x, I...)

Base.convert(::Type{GeoArray}, array::AbstractGeoArray) = GeoArray(array)

# Manually add broadcast style to GeoArray until all sources are real arrays
# and we have access to type parameter A for all of them.
# Base.BroadcastStyle(::Type{<:GeoArray{T,N,D,R,A}}) where {T,N,D,R,A} = begin
    # inner_style = typeof(Base.BroadcastStyle(A))
    # return DimensionalData.DimensionalStyle{inner_style}()
# end

# Utils ########################################################################

@inline getmeta(A::AbstractGeoArray, key, fallback) = getmeta(metadata(A), key, fallback)
@inline getmeta(m::Nothing, key, fallback) = fallback
@inline getmeta(m::Union{NamedTuple,Dict}, key, fallback) = key in keys(m) ?  m[key] : fallback
@inline getmeta(m::Metadata, key, fallback) = getmeta(val(m), key, fallback)
