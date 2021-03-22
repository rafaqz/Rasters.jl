
"""
    AbstractGeoArray <: DimensionalData.AbstractDimArray

Abstract supertype for objects that wrap an array (or location of an array) 
and metadata about its contents. It may be memory ([`GeoArray`](@ref)) or
disk-backed ([`NCDarray`](@ref), [`GDALarray`](@ref), [`GRDarray`](@ref)).

`AbstractGeoArray`s inherit from [`AbstractDimArray`]($DDarraydocs)
from DimensionalData.jl. They can be indexed as regular Julia arrays or with
DimensionalData.jl [`Dimension`]($DDdimdocs)s. They will plot as a heatmap in
Plots.jl with correct coordinates and labels, even after slicing with
`getindex` or `view`. `getindex` on a `AbstractGeoArray` will always return
a memory-backed `GeoArray`.
"""
abstract type AbstractGeoArray{T,N,D,A} <: AbstractDimensionalArray{T,N,D,A} end

# Marker singleton for lazy loaded arrays, only used for broadcasting.
# Can be removed when DiskArrays.jl is used everywhere
struct LazyArray{T,N} <: AbstractArray{T,N} end

# Interface methods ###########################################################

"""
    missingval(x)

Returns the value representing missing data in the dataset
"""
function missingval end
missingval(x) = missing
missingval(A::AbstractGeoArray) = A.missingval

"""
    crs(x)

Get the projected coordinate reference system of a `Y` or `X` `Dimension`,
or of the `Y`/`X` dims of an `AbstractGeoArray`.

For [`Mapped`](@ref) mode this may be `nothing` as there may be not projected
coordinate reference system at all.
"""
function crs end
function crs(A::AbstractGeoArray)
    if hasdim(A, Y)
        crs(dims(A, Y))
    elseif hasdim(A, X)
        crs(dims(A, X))
    else
        error("No Y or X dimension, crs not available")
    end
end
crs(dim::Dimension) = crs(mode(dim))

"""
    mappedcrs(x)

Get the mapped coordinate reference system for the `Y`/`X` dims of an array.

In [`Projected`](@ref) mode this is used to convert [`Selector`]($DDselectordocs)
values form the mappedcrs defined projection to the underlying projection, and to
show plot axes in the mapped projection.

In `Mapped` mode this is the coordinate reference system of the index values.
"""
function mappedcrs end
function mappedcrs(A::AbstractGeoArray)
    if hasdim(A, Y)
        mappedcrs(dims(A, Y))
    elseif hasdim(A, X)
        mappedcrs(dims(A, X))
    else
        error("No Y or X dimension, mappedcrs not available")
    end
end
mappedcrs(dim::Dimension) = mappedcrs(mode(dim))

# DimensionalData methods

DD.units(A::AbstractGeoArray) = getmeta(A, :units, nothing)

for f in (:mappedbounds, :projectedbounds, :mappedindex, :projectedindex)
    @eval ($f)(A::AbstractGeoArray, dims_) = ($f)(dims(A, dims_))
    @eval ($f)(A::AbstractGeoArray) = ($f)(dims(A))
end

# Rebuild all types of AbstractGeoArray as GeoArray
function DD.rebuild(
    A::AbstractGeoArray, data, dims::Tuple, refdims, name,
    metadata, missingval=missingval(A)
)
    GeoArray(data, dims, refdims, name, metadata, missingval)
end
function DD.rebuild(A::AbstractGeoArray;
    data=data(A), dims=dims(A), refdims=refdims(A), name=name(A),
    metadata=metadata(A), missingval=missingval(A)
)
    GeoArray(data, dims, refdims, name, metadata, missingval)
end

Base.parent(A::AbstractGeoArray) = data(A)


"""
Abstract supertype for all memory-backed GeoArrays where the data is an array.
"""
abstract type MemGeoArray{T,N,D,A} <: AbstractGeoArray{T,N,D,A} end

"""
    DiskGeoArray <: AbstractGeoArray

Abstract supertype for all disk-backed GeoArrays.
For these the data is lazyily loaded from disk.

To load a `DiskGeoArray` and operate on the data multiple times, use
[`open`](@ref) and a `do` block.
"""
abstract type DiskGeoArray{T,N,D,A} <: AbstractGeoArray{T,N,D,A} end

filename(A::DiskGeoArray) = A.filename

"""
    data(f, A::DiskGeoArray)

Run method `f` on the data source object for `A`, as passed by the
`withdata` method for the array. The only requirement of the
object is that it has an `Array` method that returns the data as an array.
"""
DD.data(A::DiskGeoArray) = withsourcedata(Array, A)

# Base methods

Base.size(A::DiskGeoArray) = A.size

@propagate_inbounds function Base.getindex(
    A::DiskGeoArray, i1::DD.StandardIndices, i2::DD.StandardIndices, I::DD.StandardIndices...
)
    _rebuildgetindex(A, i1, i2, I...)
end
@propagate_inbounds function Base.getindex(A::DiskGeoArray, i1::Integer, i2::Integer, I::Vararg{<:Integer})
    _rawgetindex(A, i1, i2, I...)
end
# Linear indexing returns Array
@propagate_inbounds function Base.getindex(
    A::DiskGeoArray, i::Union{Colon,AbstractVector{<:Integer}}
)
    _rawgetindex(A, i)
end
# Except 1D DimArrays
@propagate_inbounds function Base.getindex(
    A::DiskGeoArray{<:Any,1}, i::Union{Colon,AbstractVector{<:Integer}}
)
    _rebuildgetindex(A, i)
end

@propagate_inbounds function _rawgetindex(A, I...)
    withsourcedata(A) do data
        readwindowed(data, I...)
    end
end

@propagate_inbounds function _rebuildgetindex(A, I...)
    withsourcedata(A) do data
        dims_, refdims_ = DD.slicedims(dims(A), refdims(A), I)
        data = readwindowed(data, I...)
        rebuild(A, data, dims_, refdims_)
    end
end

Base.write(A::T) where T <: DiskGeoArray = write(filename(A), A)
Base.write(filename::AbstractString, A::T) where T <: DiskGeoArray = write(filename, T, A)


# Concrete implementation ######################################################

"""
    GeoArray <: MemGeoArray

    GeoArray(A::AbstractArray{T,N}, dims::Tuple; kw...)
    GeoArray(A::AbstractArray{T,N}; dims, kw...)
    GeoArray(A::AbstractGeoArray; kw...) =

A generic, memory-backed spatial array type. All [`AbstractGeoArray`](@ref) are
converted to `GeoArray` when indexed or otherwise transformed.

# Keywords

- `data`: can replace the data in an `AbstractGeoArray`
- `dims`: `Tuple` of `Dimension`s for the array.
- `refdims`: `Tuple of` position `Dimension`s the array was sliced from,
    defaulting to `()`.
- `name`: `Symbol` name for the array.
- `missingval`: Value reprsenting missing values, defaulting to `missing`.
    can be passed it.
- `metadata`: `ArrayMetadata` object for the array, or `NoMetadata()`.
"""
struct GeoArray{T,N,D<:Tuple,R<:Tuple,A<:AbstractArray{T,N},Na,Me,Mi} <: MemGeoArray{T,N,D,A}
    data::A
    dims::D
    refdims::R
    name::Na
    metadata::Me
    missingval::Mi
end
function GeoArray(A::AbstractArray, dims::Tuple;
    refdims=(), name=Symbol(""), metadata=NoMetadata(), missingval=missing
)
    GeoArray(A, DD.formatdims(A, dims), refdims, name, metadata, missingval)
end
function GeoArray(A::AbstractArray;
    dims, refdims=(), name=Symbol(""), metadata=NoMetadata(), missingval=missing
)
    GeoArray(A, DD.formatdims(A, dims), refdims, name, metadata, missingval)
end
function GeoArray(A::AbstractGeoArray;
    data=data(A), dims=dims(A), refdims=refdims(A),
    name=name(A), metadata=metadata(A), missingval=missingval(A)
)
    GeoArray(data, dims, refdims, name, metadata, missingval)
end

@propagate_inbounds function Base.setindex!(A::GeoArray, x, I::DD.StandardIndices...)
    setindex!(data(A), x, I...)
end

Base.convert(::Type{GeoArray}, array::AbstractGeoArray) = GeoArray(array)
