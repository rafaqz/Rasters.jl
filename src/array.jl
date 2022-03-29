
"""
    AbstractRaster <: DimensionalData.AbstractDimArray

Abstract supertype for objects that wrap an array (or location of an array) 
and metadata about its contents. It may be memory or hold a `FileArray`, which
holds the filename, and is only opened when required.

`AbstractRaster`s inherit from [`AbstractDimArray`]($DDarraydocs)
from DimensionalData.jl. They can be indexed as regular Julia arrays or with
DimensionalData.jl [`Dimension`]($DDdimdocs)s. They will plot as a heatmap in
Plots.jl with correct coordinates and labels, even after slicing with
`getindex` or `view`. `getindex` on a `AbstractRaster` will always return
a memory-backed `Raster`.
"""
abstract type AbstractRaster{T,N,D,A} <: AbstractDimArray{T,N,D,A} end

# Interface methods ###########################################################
"""
    missingval(x)

Returns the value representing missing data in the dataset
"""
function missingval end
missingval(x) = missing
missingval(A::AbstractRaster) = A.missingval

filename(A::AbstractRaster) = filename(parent(A))
filename(A::AbstractArray) = nothing # Fallback

cleanreturn(A::AbstractRaster) = modify(cleanreturn, A)
cleanreturn(x) = x

isdisk(A::AbstractRaster) = parent(A) isa DiskArrays.AbstractDiskArray
isdisk(x) = false
ismem(x) = !isdisk(x)

setcrs(x::AbstractRaster, crs) = set(x, setcrs(dims(x), crs)...)
setmappedcrs(x::AbstractRaster, mappedcrs) = set(x, setmappedcrs(dims(x), mappedcrs)...)

function Base.:(==)(A::AbstractRaster{T,N}, B::AbstractRaster{T,N}) where {T,N} 
    size(A) == size(B) && all(A .== B)
end

"""
    crs(x)

Get the projected coordinate reference system of a `Y` or `X` `Dimension`,
or of the `Y`/`X` dims of an `AbstractRaster`.

For [`Mapped`](@ref) lookup this may be `nothing` as there may be no projected
coordinate reference system at all.
"""
function crs end
function crs(obj)
    if hasdim(obj, Y)
        crs(dims(obj, Y))
    elseif hasdim(obj, X)
        crs(dims(obj, X))
    else
        nothing
    end
end
crs(dim::Dimension) = crs(lookup(dim))

"""
    mappedcrs(x)

Get the mapped coordinate reference system for the `Y`/`X` dims of an array.

In [`Projected`](@ref) lookup this is used to convert [`Selector`]($DDselectordocs)
values form the mappedcrs defined projection to the underlying projection, and to
show plot axes in the mapped projection.

In `Mapped` lookup this is the coordinate reference system of the index values.
"""
function mappedcrs end
function mappedcrs(obj)
    if hasdim(obj, Y)
        mappedcrs(dims(obj, Y))
    elseif hasdim(obj, X)
        mappedcrs(dims(obj, X))
    else
        nothing
    end
end
mappedcrs(dim::Dimension) = mappedcrs(lookup(dim))

for f in (:mappedbounds, :projectedbounds, :mappedindex, :projectedindex)
    @eval ($f)(A::AbstractRaster, dims_) = ($f)(dims(A, dims_))
    @eval ($f)(A::AbstractRaster) = ($f)(dims(A))
end

# DimensionalData methods

DD.units(A::AbstractRaster) = get(metadata(A), :units, nothing)
# Rebuild all types of AbstractRaster as Raster
function DD.rebuild(
    A::AbstractRaster, data, dims::Tuple, refdims, name,
    metadata, missingval=missingval(A)
)
    Raster(data, dims, refdims, name, metadata, missingval)
end
function DD.rebuild(A::AbstractRaster;
    data=parent(A), dims=dims(A), refdims=refdims(A), name=name(A),
    metadata=metadata(A), missingval=missingval(A)
)
    rebuild(A, data, dims, refdims, name, metadata, missingval)
end

function DD.modify(f, A::AbstractRaster)
    newdata = open(A) do O
        f(parent(O))
    end
    size(newdata) == size(A) || error("$f returns an array with size $(size(newdata)) when the original size was $(size(A))")
    return rebuild(A, newdata)
end

function DD.DimTable(As::Tuple{<:AbstractRaster,Vararg{<:AbstractRaster}}...)
    DimTable(DimStack(map(read, As...)))
end

# DiskArrays methods

DiskArrays.eachchunk(A::AbstractRaster) = DiskArrays.eachchunk(parent(A))
DiskArrays.haschunks(A::AbstractRaster) = DiskArrays.haschunks(parent(A))
function DA.readblock!(A::AbstractRaster, dst, r::AbstractUnitRange...)
    DA.readblock!(parent(A), dst, r...)
end
function DA.writeblock!(A::AbstractRaster, src, r::AbstractUnitRange...) 
    DA.writeblock!(parent(A), src, r...)
end

# Base methods

Base.parent(A::AbstractRaster) = A.data

"""
    open(f, A::AbstractRaster; write=false)

`open` is used to open any `AbstractRaster` and do multiple operations
on it in a safe way. The `write` keyword opens the file in write lookup so that it
can be altered on disk using e.g. a broadcast.

`f` is a method that accepts a single argument - an `Raster` object
which is just an `AbstractRaster` that holds an open disk-based object.
Often it will be a `do` block:

```julia
# A is an `Raster` wrapping the opened disk-based object.
open(Raster(filepath); write=true) do A
    mask!(A; to=maskfile)
    A[I...] .*= 2
    # ...  other things you need to do with the open file
end
```

By using a do block to open files we ensure they are always closed again
after we finish working with them.
"""
function Base.open(f::Function, A::AbstractRaster; kw...)
    # Open FileArray to expose the actual dataset object, even inside nested wrappers
    select = FileArray
    ignore = Union{Dict,Set,Base.MultiplicativeInverses.SignedMultiplicativeInverse}
    fa = Flatten.flatten(parent(A), select, ignore)
    if fa == ()
        f(Raster(parent(A), dims(A), refdims(A), name(A), metadata(A), missingval(A)))
    else
        open(fa[1]; kw...) do x
            # Rewrap the opened object where the FileArray was
            d = Flatten.reconstruct(parent(A), (x,), select, ignore) 
            f(Raster(d, dims(A), refdims(A), name(A), metadata(A), missingval(A)))
        end
    end
end
Base.write(A::T) where T <: AbstractRaster = write(filename(A), A)

# Concrete implementation ######################################################

"""
    Raster <: AbsractRaster

    Raster(filepath::AbstractString, dims; kw...)
    Raster(A::AbstractArray{T,N}, dims; kw...)
    Raster(A::AbstractRaster; kw...)

A generic [`AbstractRaster`](@ref) for spatial/raster array data. It may hold
memory-backed arrays or [`FileArray`](@ref), that simply holds the `String` path
to an unopened file. This will only be opened lazily when it is indexed with `getindex`
or when `read(A)` is called. Broadcasting, taking a view, reversing and most other 
methods _do not_ load data from disk: they are applied later, lazily.

# Keywords

- `data`: can replace the data in an `AbstractRaster`
- `dims`: `Tuple` of `Dimension`s for the array.
- `refdims`: `Tuple of` position `Dimension`s the array was sliced from, defaulting to `()`.
- `key`: `Symbol` key to desired layer in a multi-layer dataset, when a `filpath` is used.
- `name`: `Symbol` name for the array. `key` is used by default when a filepath `String` is pased in.
- `missingval`: value reprsenting missing data, normally detected form the file. Set manually
    when you know the value is not specified or is incorrect. This will *not* change any
    values in the raster, it simply assigns which value is treated as missing. To replace all of
    the missing values in the raster, use [`replace_missing`](@ref).
- `metadata`: `ArrayMetadata` object for the array, or `NoMetadata()`.
- `crs`: the coordinate reference system of  the objects `XDim`/`YDim` dimensions. 
    Only set this if you know the detected crs is incrorrect, or it is not present in
    the file.
- `mappedcrs`: the mapped coordinate reference system of the objects `XDim`/`YDim` dimensions.
    for `Mapped` lookups these are the actual values of the index. For `Projected` lookups
    this can be used to index in eg. `EPSG(4326)` lat/lon values, having it converted automatically.
    Only set this if the detected `mappedcrs` in incorrect, or the file does not have a `mappedcrs`,
    e.g. a tiff.
"""
struct Raster{T,N,D<:Tuple,R<:Tuple,A<:AbstractArray{T,N},Na,Me,Mi} <: AbstractRaster{T,N,D,A}
    data::A
    dims::D
    refdims::R
    name::Na
    metadata::Me
    missingval::Mi
end
function Raster(A::AbstractArray, dims::Tuple;
    refdims=(), name=Symbol(""), metadata=NoMetadata(), missingval=missing,
    crs=nothing, mappedcrs=nothing
)
    A = Raster(A, Dimensions.format(dims, A), refdims, name, metadata, missingval)
    A = isnothing(crs) ? A : setmappedcrs(A, crs)
    A = isnothing(mappedcrs) ? A : setmappedcrs(A, mappedcrs)
    return A
end
function Raster(A::AbstractArray{<:Any,1}, dims::Tuple{<:Dimension,<:Dimension,Vararg}; kw...)
    Raster(reshape(A, map(length, dims)), dims; kw...)
end
function Raster(table, dims::Tuple; name=first(_not_a_dimcol(table, dims)), kw...)
    Tables.istable(table) || throw(ArgumentError("First argument to `Raster` is not a table or other known object: $table"))
    isnothing(name) && throw(UndefKeywordError(:name))
    cols = Tables.columns(table)
    A = reshape(cols[name], map(length, dims))
    return Raster(A, dims; name, kw...)
end
Raster(A::AbstractArray; dims, kw...) = Raster(A, dims; kw...)
function Raster(A::AbstractDimArray;
    data=parent(A), dims=dims(A), refdims=refdims(A),
    name=name(A), metadata=metadata(A), missingval=missingval(A), kw...
)
    return Raster(data, dims; refdims, name, metadata, missingval, kw...)
end
function Raster(filename::AbstractString; name=nothing, key=name, kw...)
    _open(filename) do ds
        key = filekey(ds, key)
        Raster(ds, filename, key; kw...)
    end
end
function Raster(ds, filename::AbstractString, key=nothing;
    crs=nothing, mappedcrs=nothing, dims=nothing, refdims=(),
    name=Symbol(key isa Nothing ? "" : string(key)),
    metadata=metadata(ds), missingval=missingval(ds), write=false,
    source=_sourcetype(filename)
)
    crs = defaultcrs(source, crs)
    mappedcrs = defaultmappedcrs(source, mappedcrs)
    dims = dims isa Nothing ? DD.dims(ds, crs, mappedcrs) : dims
    data = FileArray(ds, filename; key, write)
    return Raster(data, dims, refdims, name, metadata, missingval)
end

@propagate_inbounds function Base.setindex!(A::Raster, x, I::DD.StandardIndices...)
    setindex!(parent(A), x, I...)
end

filekey(ds, key) = key
filekey(filename::String) = Symbol(splitext(basename(filename))[1])

# Precompile
precompile(Raster, (String,))

@deprecate geoarray(args...; kw...) Raster(args...; kw...)

DD.dimconstructor(::Tuple{<:Dimension{<:AbstractProjected},Vararg{<:Dimension}}) = Raster
