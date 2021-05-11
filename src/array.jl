
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

For [`Mapped`](@ref) mode this may be `nothing` as there may be no projected
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
    rebuild(A, data, dims, refdims, name, metadata, missingval)
end


Base.parent(A::AbstractGeoArray) = data(A)

filename(A::AbstractGeoArray) = filename(data(A))

Base.write(A::T) where T <: AbstractGeoArray = write(filename(A), A)

# Concrete implementation ######################################################

"""
    GeoArray <: AbsractGeoArray

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
struct GeoArray{T,N,D<:Tuple,R<:Tuple,A<:AbstractArray{T,N},Na,Me,Mi} <: AbstractGeoArray{T,N,D,A}
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
    data=readwindowed(data(A)), dims=dims(A), refdims=refdims(A),
    name=name(A), metadata=metadata(A), missingval=missingval(A)
)
    GeoArray(data, dims, refdims, name, metadata, missingval)
end
function GeoArray(filename::AbstractString; key=nothing, kw...)
    isfile(filename) || error("File not found: $filename")
    _read(filename) do ds
        GeoArray(ds, filename, key; kw...)
    end
end
function GeoArray(ds, filename::AbstractString, key=nothing;
    crs=nothing, mappedcrs=nothing,
    dims=dims(ds, crs, mappedcrs), refdims=(),
    name=Symbol(key isa Nothing ? "" : string(key)),
    metadata=metadata(ds),
    missingval=missingval(ds),
)
    source = _sourcetype(filename)
    crs = defaultcrs(source, crs)
    mappedcrs = defaultmappedcrs(source, mappedcrs)
    data = FileArray(ds, filename, key)
    GeoArray(data, dims, refdims, name, metadata, missingval)
end

@propagate_inbounds function Base.setindex!(A::GeoArray, x, I::DD.StandardIndices...)
    setindex!(data(A), x, I...)
end

filename(A::GeoArray) = filename(data(A))
