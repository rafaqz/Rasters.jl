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

# Interface methods ############################################################

for func in (:dims, :metadata, :missingval)
    @eval begin
        $func(s::AbstractGeoStack) = $func(s, first(keys(s)))
        $func(s::AbstractGeoStack, key::Key) =
            safeapply(dataset -> $func(s, dataset, key), s, source(s, key))
        $func(s::AbstractGeoStack, source, key::Key) = $func(source)
    end
end
refdims(s::AbstractGeoStack) = s.refdims
window(s::AbstractGeoStack) = s.window


"""
    source(s::AbstractGeoStack, [key])

The data `source` of an AbstractGeoStack is an abstraction that can be either:
a named tuple of filepaths to files containing single layers; a single filepath
for a file containing layers with multiple keys; or a NamedTuple that contains
multiple AbstractGeoArrays.

It allows the same methods to apply to many different file types and structures.
"""
function source end

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

Base.copy(stack::AbstractGeoStack) = rebuild(stack; data=map(copy, source(stack)))

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
Base.keys(s::AbstractGeoStack{<:AbstractString}) = Symbol.(safeapply(keys, s, source(s)))
Base.keys(s::AbstractGeoStack{<:NamedTuple}) = Symbol.(keys(source(s)))
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


# abstract disk-based stack ######################################################

"""
Abstract supertype for all disk backed [`AbstractGeoStack`](@ref)s
"""
abstract type DiskGeoStack{T} <: AbstractGeoStack{T} end

source(s::DiskGeoStack, args...) = filename(s, args...)

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

@inline Base.view(s::DiskGeoStack, I...) = rebuild(s; window=I)


# abstract memory-based stack ######################################################

"""
[`AbstractGeoStack`](@ref) backed by memory.
"""
abstract type MemGeoStack{T} <: AbstractGeoStack{T} end

data(s::MemGeoStack) = s.data
data(s::MemGeoStack{<:NamedTuple}, key::Key) = data(s)[Symbol(key)]

source(s::MemGeoStack{<:NamedTuple}, args...) = data(s, args...)

@inline Base.view(s::MemGeoStack, I...) =
    rebuild(s; data=NamedTuple{keys(s)}(view(a, I...) for a in values(s)))

# Concrete MemGeoStack implementation ######################################################

"""
Concrete `MemGeoStack` implementation. Holds concrete [`GeoArray`](@ref) layers in memory.

- `data`: A `NamedTuple` of `GeoArray`.
- `window`: A `Tuple` of Dimensions/Selectors/Indices that will be applied to the contained
            arrays when they are accessed.
- `refdims`: reference dimensions from earlier subsetting.
- `metadata`: Any metadata object.
"""
struct GeoStack{T,R,W,M} <: MemGeoStack{T}
    data::T
    refdims::R
    window::W
    metadata::M
end
"""
    GeoStack(data::Vararg{<:AbstractGeoArray}; kwargs...)

Convert `GeoArray`s arguments to a `GeoStack`.
"""
GeoStack(data::AbstractGeoArray...; keys=name.(data), kwargs...) =
    GeoStack(NamedTuple{_stackkeys(keys)}(data); kwargs...)
"""
    GeoStack(data::NamedTuple; [refdims=()], [window=()], [metadata=nothing]) =

Construct a `GeoStack` from a NamedTuple of [`GeoArray`](@ref) and keyword arguments.
"""
GeoStack(data::NamedTuple; refdims=(), window=(), metadata=nothing) =
    GeoStack(data, refdims, window, metadata)
"""
    GeoStack(s::AbstractGeoStack; [keys, data, refdims, window, metadata])

Construct a `GeoStack` from another `GeoStack` and keyword arguments.
`data` is a `NamedTuple` of `GeoArray`.
"""
GeoStack(s::AbstractGeoStack; keys=_stackkeys(Base.keys(s)),
         data=NamedTuple{keys}((GeoArray(s[key]) for key in keys)),
         refdims=refdims(s), window=window(s), metadata=metadata(s)) =
    GeoStack(data, refdims, window, metadata)

_stackkeys(keys) = Tuple(Symbol.(keys))

metadata(s::GeoStack) = s.metadata

safeapply(f, ::GeoStack, source) = f(source)

@inline Base.getindex(s::GeoStack, key::Key, i1::StandardIndices, I::StandardIndices...) =
    data(s)[key][i1, I...]
@inline Base.getindex(s::GeoStack, key::Key) = data(s)[key]

Base.convert(::Type{GeoStack}, src::AbstractGeoStack) = GeoStack(src)
