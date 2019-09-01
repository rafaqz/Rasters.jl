
const Key = Union{Symbol,AbstractString}

"""
Stack object for holding multiple arrays and datasets
with the same spatial metadata and bounds. As in R's raster stack.

Contained objects should share some common dims D?
"""
abstract type AbstractGeoStack{T} end

@inline data(stack, key, I...) =
    run(ds -> data(stack, ds, key, I...), stack, source(stack, key))

# We use `do` blocks to get dims and metadata from the source file safely
for func in (:dims, :metadata, :missingval)
    @eval begin
        @inline $func(stack::AbstractGeoStack) = $func(stack, first(keys(stack)))
        @inline $func(stack::AbstractGeoStack{<:NamedTuple}, key::Key) =
            run(ds -> $func(stack, ds, key), stack, source(stack, key))
        @inline $func(stack::AbstractGeoStack{<:AbstractString}, key::Key) =
            run(ds -> $func(stack, ds, key), stack, source(stack, key))
        @inline $func(stack::AbstractGeoStack, source, key::Key) = $func(source)
    end
end
@inline refdims(s::AbstractGeoStack) = s.refdims
@inline metadata(s::AbstractGeoStack) = s.metadata

# Named tuple of paths to single-layer files
@inline source(stack::AbstractGeoStack{<:NamedTuple}) = parent(stack)
@inline source(stack::AbstractGeoStack{<:NamedTuple}, key::Key) = parent(stack)[Symbol(key)]
# Single filepath for multi-layer files
@inline source(stack::AbstractGeoStack{<:AbstractString}, args...) = parent(stack)



@inline rebuild(s::AbstractGeoStack, values::Base.Generator) =
    rebuild(s, NamedTuple{keys(s)}(values))
@inline rebuild(s::AbstractGeoStack, data, newdims, newref=refdims(s)) =
    GeoStack(data, newdims, newref, metadata(s))

Base.parent(s::AbstractGeoStack) = s.data

Base.copy!(dst::AbstractGeoStack, src::AbstractGeoStack, destkeys=keys(dst)) = begin
    for key in destkeys
        key in Symbol.(keys(dst)) || throw(ArgumentError("key $key not found in dest keys"))
        key in Symbol.(keys(src)) || throw(ArgumentError("key $key not found in source keys"))
    end
    for key in Symbol.(keys(dst))
        copy!(dst[key], src, key)
    end
end

# Array interface
@inline Base.view(s::AbstractGeoStack, I...) = rebuild(s, (view(a, I...) for a in values(s)))
@inline Base.getindex(s::AbstractGeoStack, I...) = rebuild(s, (a[I...] for a in values(s)))
# Dict/Array hybrid
@inline Base.getindex(stack::AbstractGeoStack, key::Key, I::Vararg{<:AbstractDimension}) =
    getindex(stack, key, dims2indices(dims(stack, key), I)...)
@inline Base.getindex(stack::AbstractGeoStack, key::Key, I...) = data(stack, key, I...)
@inline Base.getindex(stack::AbstractGeoStack, key::Key) = data(stack, key)

@inline Base.values(stack::AbstractGeoStack) = (stack[key] for key in keys(stack))
@inline Base.length(stack::AbstractGeoStack) = length(keys(stack))
@inline Base.keys(stack::AbstractGeoStack{<:AbstractString}) = Symbol.(run(keys, stack, source(stack)))
@inline Base.keys(stack::AbstractGeoStack) = Symbol.(keys(parent(stack)))
@inline Base.names(stack::AbstractGeoStack) = keys(stack)

@mix struct GeoStackMixin{T,D,R,M}
    data::T
    dims::D
    refdims::R
    metadata::M
end

"""
Basic stack object. Holds AbstractGeoArray layers

`view` or `getindex` return another stack with the method
applied to all layers.
"""
@GeoStackMixin struct GeoStack{T,D,R,M} <: AbstractGeoStack{T} end

GeoStack(data::Vararg{<:AbstractGeoArray};
         keys=Symbol.(name.(data)), dims=(), refdims=(first(data)), metadata=nothing) =
    GeoStack(NamedTuple{keys}(data), dims, refdims, metadata)
GeoStack(stack::AbstractGeoStack; keys=keys(stack), dims=(), 
         refdims=refdims(stack), metadata=metadata(stack)) = begin
    keys = Tuple(Symbol.(keys))
    data = NamedTuple{keys}((GeoArray(stack[key]) for key in keys))
    # GeoStack(data, dims, refdims, metadata)
end

# GeoStack is in-memory so we don't have to fetch anything here,
# just run the function on the contained object.
@inline run(f, stack::GeoStack, source) = f(source)
@inline data(stack::GeoStack, source, key::Key, I...) = source[I...]

# GeoStack keys are in memory objects so we just return them
@inline Base.getindex(stack::GeoStack, key::Key) = parent(stack)[Symbol(key)]
@inline Base.getindex(stack::GeoStack, key::Key) = getindex(parent(stack), key)

Base.convert(::GeoStack, stack::AbstractGeoStack) = GeoStack(stack)
