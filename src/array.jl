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

missingval(A::AbstractGeoArray) = A.missingval

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

# DimensionalData methods

units(A::AbstractGeoArray) = getmeta(A, :units, "")
label(A::AbstractGeoArray) = string(name(A), " ", units(A))

# Rebuild all types of AbstractGeoArray as GeoArray
rebuild(A::AbstractGeoArray, data, dims::Tuple, refdims, name, metadata, missingval=missingval(A)) =
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

"""
    data(f, A::DiskGeoArray)

Run method `f` on the data source object for A, as passed by the
`withdata` method for the array. The only requirement of the
object is that it has an `Array` method that returns the data as an array.
"""
data(A::DiskGeoArray) = withsourcedata(Array, A)

# Base methods

Base.size(A::DiskGeoArray) = A.size

Base.getindex(A::DiskGeoArray, I::Vararg{<:Union{<:Integer,<:AbstractArray}}) =
    withsourcedata(A) do data
        dims_, refdims_ = slicedims(dims(A), refdims(A), I)
        data = readwindowed(data, I...)
        rebuild(A, data, dims_, refdims_)
    end
Base.getindex(A::DiskGeoArray, i1::Integer, I::Vararg{<:Integer}) =
    withsourcedata(A) do data
        readwindowed(data, i1, I...)
    end

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
@inline GeoArray(A::AbstractArray, dims::Tuple;
                 refdims=(), name="", metadata=nothing, missingval=missing) =
    GeoArray(A, formatdims(A, dims), refdims, name, metadata, missingval)
"""
    GeoArray(A::AbstractGeoArray; [data=data(A), dims=dims(A), refdims=refdims(A),
             name=name(A), metadata=metadata(A), missingval=missingval(A)]) =

Construct a [`GeoArray`](@ref) from another [`AbstractGeoArray`](@ref), and
keyword arguments.
"""
@inline GeoArray(A::AbstractGeoArray; data=data(A), dims=dims(A), refdims=refdims(A),
                 name=name(A), metadata=metadata(A), missingval=missingval(A)) =
    GeoArray(data, dims, refdims, name, metadata, missingval)

Base.@propagate_inbounds Base.setindex!(A::GeoArray, x, I::Vararg{DimensionalData.StandardIndices}) =
    setindex!(data(A), x, I...)

Base.convert(::Type{GeoArray}, array::AbstractGeoArray) = GeoArray(array)
