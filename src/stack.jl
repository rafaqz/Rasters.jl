
const Key = Union{Symbol,AbstractString}

"""
Stack object for holding multiple arrays and datasets
with the same spatial metadata and bounds. As in R's raster stack.

Contained objects should share some common dims D?
"""
abstract type AbstractGeoStack{T} end

(::Type{T})(data, keys; kwargs...) where T<:AbstractGeoStack =
    T(NamedTuple{Tuple(Symbol.(keys))}(Tuple(data)); kwargs...)
(::Type{T})(data; refdims=(), window=()) where T<:AbstractGeoStack =
    T(data, window, refdims)

@premix struct GeoStackMixin{T,W,R}
    data::T
    window::W
    refdims::R
end

# Get the data source with a wrapper method. This means disk datasets are
# correctly closed after use but are accessed in a standardised way.
data(s, key, I...) =
    run(ds -> data(s, ds, key, I...), s, source(s, key))

# Core methods
for func in (:dims, :metadata, :missingval)
    @eval begin
        $func(s::AbstractGeoStack) = $func(s, first(keys(s)))
        $func(s::AbstractGeoStack, key::Key) =
            run(ds -> $func(s, ds, key), s, source(s, key))
        $func(s::AbstractGeoStack, source, key::Key) = $func(source)
    end
end
window(s::AbstractGeoStack) = s.window
refdims(s::AbstractGeoStack) = s.refdims

# The data source is abstracted to handle a named tuple of filepaths,
# a single filepath containing layers with multiple keys, or a NamedTuple
# that contains memory back GeoArrays
source(s::AbstractGeoStack{<:NamedTuple}) = parent(s)
source(s::AbstractGeoStack{<:NamedTuple}, key::Key) = parent(s)[Symbol(key)]
# Single filepath for multi-layer files
source(s::AbstractGeoStack{<:AbstractString}, args...) = parent(s)



rebuild(s::AbstractGeoStack, values::Base.Generator; kwargs...) =
    rebuild(s; parent=NamedTuple{keys(s)}(values), kwargs...)
rebuild(s::AbstractGeoStack; parent=parent(s), refdims=refdims(s),
        window=window(s), metadata=metadata(s)) =
    GeoStack(parent, refdims, window, metadata)

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
Base.copy!(dst::AbstractArray, src::AbstractGeoStack, key) = copy!(dst, src[key])

# Array interface: delete this and just use window?
@inline Base.view(s::AbstractGeoStack, I...) = rebuild(s, (view(a, I...) for a in values(s)))
@inline Base.getindex(s::AbstractGeoStack, I...) = rebuild(s, (a[I...] for a in values(s)))

# Dict/Array hybrid
Base.getindex(s::AbstractGeoStack, key::Key, I::Vararg{<:AbstractDimension}) =
    getindex(s, key, dims2indices(dims(s, key), I)...)
Base.getindex(s::AbstractGeoStack, key::Key, I...) = 
    data(s, key, applywindow(s, key, I)...)
Base.getindex(s::AbstractGeoStack, key::Key) = data(s, key, windoworempty(s, key)...)

Base.values(s::AbstractGeoStack) = (s[key] for key in keys(s))
Base.length(s::AbstractGeoStack) = length(keys(s))
Base.keys(s::AbstractGeoStack{<:AbstractString}) = Symbol.(run(keys, s, source(s)))
Base.keys(s::AbstractGeoStack) = Symbol.(keys(parent(s)))
Base.names(s::AbstractGeoStack) = keys(s)


"""
Basic stack object. Holds AbstractGeoArray layers

`view` or `getindex` return another stack with the method
applied to all layers.
"""
@GeoStackMixin struct GeoStack{M} <: AbstractGeoStack{T}
    metadata::M
end

stackkeys(keys) = Tuple(Symbol.(keys))

GeoStack(data::Vararg{<:AbstractGeoArray}; keys=name.(data), kwargs...) = 
    GeoStack(NamedTuple{stackkeys(keys)}(data); kwargs...)
GeoStack(data::NamedTuple; refdims=(), window=(), metadata=nothing) =
    GeoStack(data, window, refdims, metadata)
GeoStack(s::AbstractGeoStack;
         keys=stackkeys(Base.keys(s)),
         parent = NamedTuple{keys}((GeoArray(s[key]) for key in keys)),
         dims=dims(s), refdims=refdims(s),
         metadata=metadata(s), window=window(s)) =
    GeoStack(parent, window, refdims, metadata)

metadata(s::GeoStack) = s.metadata
rebuild(s::GeoStack; parent=parent(s), refdims=refdims(s),
        window=window(s), metadata=metadata(s)) =
    GeoStack(parent, window, refdims, metadata)

# GeoStack is in-memory so we don't have to fetch anything here,
# just run the function on the contained object.
run(f, s::GeoStack, source) = f(source)
data(s::GeoStack, key::Key, I...) = data(s, key)[I...]
data(s::GeoStack, key::Key) = parent(s)[key]

# GeoStack keys are in-memory objects so we just return them
# @inline Base.getindex(s::GeoStack, key::Key) = parent(s)[Symbol(key)]
# @inline Base.getindex(s::GeoStack, key::Key) = getindex(parent(s), key)

Base.convert(::GeoStack, s::AbstractGeoStack) = GeoStack(s)
