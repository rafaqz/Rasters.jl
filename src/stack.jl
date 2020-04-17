# Accept either Symbol or String keys, but allways convert to Symbol
const Key = Union{Symbol,AbstractString}

"""
Stack objects hold multiple raster array that share spatial metadata and bounds.

These are NamedTuple-like structures that may either contain `NamedTuple`
of `AbstractGeoArray`, string paths that will load `AbstractGeoArray`, or a single
path that points to as a multi-layered stack of arrays. 

The primary purpose of  is that use and syntax is identical for all cases,
abstracting away data source and simplifying access code. `getindex` on any
`AbstractGeoStack` may return a memory backed standard `GeoArray`, or a disk
base AbstractGeoArray. `geoarray[:somelayer] |> plot` plots the layers array,
while `geoarray[:somelayer, Lon(1:100), Band(2)] |> plot` will plot the
subsetted array directly from disk, without loading the whole array. 

`getindex` on a GeoStack returns another stack with the method applied to all layers.
"""
abstract type AbstractGeoStack{T} end

(::Type{T})(data, keys; kwargs...) where T<:AbstractGeoStack =
    T(NamedTuple{Tuple(Symbol.(keys))}(Tuple(data)); kwargs...)

# Standard fields

refdims(s::AbstractGeoStack) = s.refdims
metadata(s::AbstractGeoStack) = s.metadata
window(s::AbstractGeoStack) = s.window
"""
    kwargs(s::AbstractGeoStack)

Returns the keyword arguments that will be passed to the child array constructor.
"""
kwargs(s::AbstractGeoStack) = s.kwargs

# Interface methods ############################################################

# Default dims, metadata and missingval methods
#
# For DiskStack we query the underlying object - avoiding building
# an AbstractGeoArray unless we have to. Examples are be an NCDatasets,
# ArchGDAL, or HDF5 `Dataset` object. These sources will add a 
# For MemGeoStack the childobj is just a GeoArray as it is allready loaded in memory.
for func in (:dims, :metadata, :missingval)
    @eval begin
        $func(s::AbstractGeoStack) = $func(s, first(keys(s)))
        $func(s::AbstractGeoStack, key::Key) =
            querychild(childobj -> $func(s, childobj, key), childsource(s, key), s)
        $func(s::AbstractGeoStack, childobj, key::Key) = $func(childobj)
    end
end


"""
    childsourc(s::AbstractGeoStack, [key])

Get the lower lovel child object. This can be an `AbstractGeoArray` or
a lower-level object with GeoData methods defined. Returning
the low-level object can be better performance as we do not have to
processes everything needed to build a full `AbstractGeoArray`.
"""
function childsource end

# Base methods #################################################################

"""
    Base.write(filename::AbstractString, T::Type{<:AbstractGeoArray}, s::AbstractGeoStack)

Save all layers of an `AbstractGeoStack` to separate files, using the backend determined
by `T`.

## Example

```julia
write(filename, GDALArray, stack)
```
"""
Base.write(filename::AbstractString, ::Type{T}, s::AbstractGeoStack) where T <: AbstractGeoArray =
    for key in keys(s)
        base, ext = splitext(filename)
        fn = joinpath(string(base, "_", key, ext))
        write(fn, T, s[key])
    end

Base.copy(stack::AbstractGeoStack) = rebuild(stack; data=map(copy, childsource(stack)))

"""
    Base.copy!(dst::AbstractGeoStack, src::AbstractGeoStack, [keys=keys(dst)])

Copy all or a subset of layers from one stack to another.

## Example
```julia
copy!(dst::AbstractGeoStack, src::AbstractGeoStack, keys=(:sea_surface_temp, :humidity))
```
"""
Base.copy!(dst::AbstractGeoStack, src::AbstractGeoStack, keys=keys(dst)) = begin
    for key in keys
        key in Symbol.(Base.keys(dst)) || throw(ArgumentError("key $key not found in dest keys"))
        key in Symbol.(Base.keys(src)) || throw(ArgumentError("key $key not found in source keys"))
    end
    for key in Symbol.(keys)
        copy!(dst[key], src, key)
    end
end
Base.copy!(dst::AbstractArray, src::AbstractGeoStack, key) =
    copy!(dst, src[key])

@inline Base.getindex(s::AbstractGeoStack, I...) =
    rebuild(s; data=NamedTuple{keys(s)}(a[I...] for a in values(s)))

# Dict/Array hybrid
@inline Base.getindex(s::AbstractGeoStack, key::Key, I::Vararg{<:Dimension}) =
    getindex(s, key, dims2indices(dims(s), I)...)

Base.values(s::AbstractGeoStack) = (s[key] for key in keys(s))
Base.length(s::AbstractGeoStack) = length(keys(s))
Base.keys(s::AbstractGeoStack{<:AbstractString}) = Symbol.(querychild(keys, childsource(s), s))
Base.keys(s::AbstractGeoStack{<:NamedTuple}) = Symbol.(keys(childsource(s)))
Base.names(s::AbstractGeoStack) = keys(s)

"""
    Base.cat(stacks::AbstractGeoStack...; [keys=keys(stacks[1])], dims)

Concatenate all or a subset of layers for all passed in stacks.

## Example
```julia
cat(stacks...; keys=(:sea_surface_temp, :humidity), dims)
```
"""
Base.cat(stacks::AbstractGeoStack...; keys=keys(stacks[1]), dims) = begin
    vals = Tuple(cat((s[key] for s in stacks)...; dims=dims) for key in keys)
    GeoStack(stacks[1], data=NamedTuple{keys}(vals))
end
Base.first(s::AbstractGeoStack) = s[first(keys(s))]
Base.last(s::AbstractGeoStack) = s[last(keys(s))]


# Memory-based stacks ######################################################

"""
[`AbstractGeoStack`](@ref) backed by memory.
"""
abstract type MemGeoStack{T} <: AbstractGeoStack{T} end

data(s::MemGeoStack) = s.data
data(s::MemGeoStack{<:NamedTuple}, key::Key) = data(s)[Symbol(key)]

childsource(s::MemGeoStack{<:NamedTuple}, args...) = data(s, args...)

querychild(f, childobj, ::AbstractGeoStack) = f(childobj)

childdata(f, childobj, ::AbstractGeoStack) = f(childobj)

@inline Base.view(s::MemGeoStack, I...) =
    rebuild(s; data=NamedTuple{keys(s)}(view(a, I...) for a in values(s)))

@inline Base.getindex(s::MemGeoStack, key::Key) = 
    rebuild(data(s)[key]; window=window(s), kwargs(s)...)
@inline Base.getindex(s::MemGeoStack, key::Key, i1::StandardIndices, I::StandardIndices...) =
    rebuild(data(s)[key]; window=window(s), kwargs(s)...)[i1, I...]


# Concrete MemGeoStack implementation ######################################################

"""
Concrete `MemGeoStack` implementation. Holds concrete [`GeoArray`](@ref) layers in memory.

- `data`: A `NamedTuple` of `GeoArray`.
- `window`: A `Tuple` of Dimensions/Selectors/Indices that will be applied to the contained
            arrays when they are accessed.
- `refdims`: reference dimensions from earlier subsetting.
- `metadata`: Any metadata object.
"""
struct GeoStack{T,R,W,M,K} <: MemGeoStack{T}
    data::T
    refdims::R
    window::W
    metadata::M
    kwargs::K
end
"""
    GeoStack(data::Vararg{<:AbstractGeoArray}; kwargs...)

Convert `GeoArray`s to a `GeoStack`.
"""
GeoStack(data::AbstractGeoArray...; keys=name.(data), kwargs...) =
    GeoStack(NamedTuple{cleankeys(keys)}(data); kwargs...)
"""
    GeoStack(data::NamedTuple; [window=()], [metadata=nothing], kwargs...) =

Construct a `GeoStack` from a NamedTuple of [`GeoArray`](@ref) and keyword arguments.
"""
GeoStack(data::NamedTuple; refdims=(), window=(), metadata=nothing, kwargs...) =
    GeoStack(data, refdims, window, metadata, kwargs)
"""
    GeoStack(s::AbstractGeoStack; [keys, data, refdims, window, metadata])

Construct a `GeoStack` from another `GeoStack` and keyword arguments.
`data` is a `NamedTuple` of `GeoArray`.
"""
GeoStack(s::AbstractGeoStack; keys=cleankeys(Base.keys(s)),
         data=NamedTuple{keys}((GeoArray(s[key]) for key in keys)),
         refdims=refdims(s), window=window(s), metadata=metadata(s), kwargs...) =
    GeoStack(data, refdims, window, metadata, kwargs)

Base.convert(::Type{GeoStack}, src::AbstractGeoStack) = GeoStack(src)


# Disk-based stacks ######################################################

"""
Abstract supertype for all disk backed [`AbstractGeoStack`](@ref)s
"""
abstract type DiskGeoStack{T} <: AbstractGeoStack{T} end

childsource(s::DiskGeoStack, args...) = filename(s, args...)

childtype(stack::DiskGeoStack) = stack.childtype

"""
    filename(s::DiskGeoStack)

Return the filename field of a `DiskGeoStack`. This may be a `Vector` of `String`,
or a `String`.
"""
filename(s::DiskGeoStack) = s.filename
"""
    filename(s::DiskGeoStack, key)

Return the filename field of a `DiskGeoStack` for a given key.

This will always be a single string. However, in some cases 
all keys may have the same filename.
"""
filename(s::DiskGeoStack{<:NamedTuple}, key::Key) = filename(s)[Symbol(key)]
filename(s::DiskGeoStack, key::Key) = filename(s)

# AbstractGeoStack methods

querychild(f, path::AbstractString, s::DiskGeoStack) = 
    querychild(f, childtype(s), path, s)

# Base methods

Base.getindex(s::DiskGeoStack, key::Key) =
    querychild(childtype(s), filename(s, key), s) do querychild
        childtype(s)(querychild; name=string(key), window=window(s), kwargs(s)...)
    end
# TODO make this more memory efficient - the window and index should both
# be passed to reading the array so the smallest possible window is read
# from disk
Base.getindex(s::DiskGeoStack, key::Key, i1::StandardIndices, I::StandardIndices...) = 
    s[key][i1, I...]

"""
    Base.copy!(dst::AbstractArray, src::GDALstack, key::Key)

Copy the stack layer `key` to `dst`, which can be any `AbstractArray`.
"""
Base.copy!(dst::AbstractArray, src::DiskGeoStack, key::Key) =
    querychild(childtype(src), filename(src, key), src) do dataset
        key = string(key)
        _window = maybewindow2indices(dataset, dims(dataset), window(src))
        copy!(dst, readwindowed(dataset, _window))
    end
Base.copy!(dst::AbstractGeoArray, src::DiskGeoStack, key::Key) =
    copy!(data(dst), src, key)

@inline Base.view(s::DiskGeoStack, I...) = rebuild(s; window=I)

# Concrete DiskGeoStack implementation ######################################################

"""
    DiskStack(filename::NamedTuple; window=())

Load a stack of files lazily from disk.

# Arguments
- `filename`: a NamedTuple of `String` filenames.

# Keyword arguments
- `window`: can be a tuple of Dimensions, selectors or regular indices.
- `childtype`: the type of the child data. eg. `GDALarray`.
- `kwargs`: keyword arguments to pass to the child constructor
"""
struct DiskStack{T,R,W,M,C,K} <: DiskGeoStack{T}
    filename::T
    refdims::R
    window::W
    metadata::M
    childtype::C
    kwargs::K
end
DiskStack(filename::NamedTuple; refdims=(), window=(), metadata=nothing, childtype, kwargs...) =
    DiskStack(filename, refdims, window, metadata, childtype, kwargs)


# Utils

cleankeys(keys) = Tuple(Symbol.(keys))

