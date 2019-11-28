# Accept either Symbol or String keys, but allways convert to Symbol
const Key = Union{Symbol,AbstractString}

"""
Stack objects hold multiple raster array that 
share spatial metadata and bounds. 
"""
abstract type AbstractGeoStack{T} end

(::Type{T})(data, keys; kwargs...) where T<:AbstractGeoStack =
    T(NamedTuple{Tuple(Symbol.(keys))}(Tuple(data)); kwargs...)

abstract type DiskGeoStack{T} <: AbstractGeoStack{T} end

abstract type MemGeoStack{T} <: AbstractGeoStack{T} end

# Interface methods ############################################################

for func in (:dims, :metadata, :missingval)
    @eval begin
        $func(s::AbstractGeoStack) = $func(s, first(keys(s)))
        $func(s::AbstractGeoStack, key::Key) =
            safeapply(dataset -> $func(s, dataset, key), s, source(s, key))
        $func(s::AbstractGeoStack, source, key::Key) = $func(source)
    end
end
dims(s::AbstractGeoStack) = s.dims
refdims(s::AbstractGeoStack) = s.refdims
window(s::AbstractGeoStack) = s.window

@inline rebuild(s::AbstractGeoStack, values::Base.Generator; kwargs...) =
    rebuild(s; data=NamedTuple{keys(s)}(values), kwargs...)
@inline rebuild(s::AbstractGeoStack; data=parent(s), dims=dims(s), refdims=refdims(s),
        window=window(s), metadata=metadata(s)) =
    GeoStack(data, dims, refdims, window, metadata)

"""
    data(stack, key, I...)

Get the data source with key and optionalindex using a wrapper method.

This means disk-based datasets are closed correctly after use but are
accessed in a standardised way, abstracting the actual calls for various
types.

Extending types should implement various combinations of
`data(stack::NewType, dataset, key, I...)` to efficiently load
either whole or partial datasets, depending on the arguments.
"""
data(stack, key, I...) =
    safeapply(dataset -> data(stack, dataset, key, I...), stack, source(stack, key))

"""
    source(s::AbstractGeoStack, [key])

The data `source` of an AbstractGeoStack is an abstraction that can be either: 
a named tuple of filepaths to files containing single layers; a single filepath 
for a file containing layers with multiple keys; or a NamedTuple that contains 
multiple AbstractGeoArrays.

This abstraction allows the same methods to apply to many different file types 
and structures.
"""
source(s::AbstractGeoStack{<:NamedTuple}) = parent(s)
source(s::AbstractGeoStack{<:NamedTuple}, key::Key) = parent(s)[Symbol(key)]
# Single filepath for multi-layer files
source(s::AbstractGeoStack{<:AbstractString}, args...) = parent(s)


# Base methods #################################################################

Base.parent(s::AbstractGeoStack) = s.data

Base.copy!(dst::AbstractGeoStack, src::AbstractGeoStack, destkeys=keys(dst)) = begin
    for key in destkeys
        key in Symbol.(keys(dst)) || throw(ArgumentError("key $key not found in dest keys"))
        key in Symbol.(keys(src)) || throw(ArgumentError("key $key not found in source keys"))
    end
    for key in Symbol.(destkeys)
        copy!(dst[key], src, key)
    end
end
Base.copy!(dst::AbstractArray, src::AbstractGeoStack, key) = copy!(dst, parent(src[key]))

# Array interface: delete this and just use window?
@inline Base.view(s::AbstractGeoStack, I...) = rebuild(s, (view(a, I...) for a in values(s)))
@inline Base.getindex(s::AbstractGeoStack, I...) = rebuild(s, (a[I...] for a in values(s)))

# Dict/Array hybrid
@inline Base.getindex(s::AbstractGeoStack, key::Key, I::Vararg{<:AbstractDimension}) =
    getindex(s, key, dims2indices(dims(s), I)...)
@inline Base.getindex(s::AbstractGeoStack, key::Key, I...) =
    data(s, key, applywindow(s, key, I)...)
@inline Base.getindex(s::AbstractGeoStack, key::Key) = data(s, key, windoworempty(s)...)

Base.values(s::AbstractGeoStack) = (s[key] for key in keys(s))
Base.length(s::AbstractGeoStack) = length(keys(s))
Base.keys(s::AbstractGeoStack{<:AbstractString}) = Symbol.(safeapply(keys, s, source(s)))
Base.keys(s::AbstractGeoStack) = Symbol.(keys(parent(s)))
Base.names(s::AbstractGeoStack) = keys(s)


# Concrete implementation ######################################################

"""
Basic stack object. Holds concrete GeoArray layers.

`view` or `getindex` return another stack with the method
applied to all layers.
"""
struct GeoStack{T,D,R,W,M} <: AbstractGeoStack{T}
    data::T
    dims::D
    refdims::R
    window::W
    metadata::M
end

stackkeys(keys) = Tuple(Symbol.(keys))

GeoStack(data::Vararg{<:AbstractGeoArray}; keys=name.(data), kwargs...) =
    GeoStack(NamedTuple{stackkeys(keys)}(data); kwargs...)
GeoStack(data::NamedTuple; 
         dims=dims(first(values(data))),
         refdims=refdims(first(values(data))),
         window=(), metadata=nothing) =
    GeoStack(data, dims, refdims, window, metadata)
GeoStack(s::AbstractGeoStack;
         keys=stackkeys(Base.keys(s)),
         data=NamedTuple{keys}((GeoArray(s[key]) for key in keys)),
         dims=dims(first(data)), 
         refdims=refdims(first(data)), 
         window=(),
         metadata=metadata(s)) =
    GeoStack(data, dims, refdims, window, metadata)

metadata(s::GeoStack) = s.metadata
rebuild(s::GeoStack; data=parent(s), dims=dims(s), refdims=refdims(s),
        window=window(s), metadata=metadata(s)) =
    GeoStack(data, dims, refdims, window, metadata)

# GeoStack is in-memory so we don't have to fetch anything here,
# just run the function on the contained object.
safeapply(f, ::GeoStack, data) = f(data)
data(s::GeoStack, key::Key, I...) = data(s, key)[I...]
data(s::GeoStack, key::Key) = parent(s)[key]

# GeoStack keys are in-memory objects so we just return them
# @inline Base.getindex(s::GeoStack, key::Key) = parent(s)[Symbol(key)]
# @inline Base.getindex(s::GeoStack, key::Key) = getindex(parent(s), key)

Base.convert(::Type{GeoStack}, src::AbstractGeoStack) = GeoStack(from)
