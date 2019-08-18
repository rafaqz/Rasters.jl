"""
Stack object for holding multiple arrays and datasets
with the same spatial metadata and bounds. As in Rs raster stack.

Contained objects should share common dims D?
"""
abstract type AbstractGeoStack end

dims(s::AbstractGeoStack) = s.dims
metatdata(s::AbstractGeoStack) = s.metatdata
keysource(s::AbstractGeoStack) = parent(s)

# Dict/NamedTuple interface
Base.keys(s::AbstractGeoStack) = keys(keysource(s))
Base.values(s::AbstractGeoStack) = values(keysource(s))
Base.length(s::AbstractGeoStack) = length(keysource(s))

""" 
Basic stack object. Holds AbstractGeoArray layers 

`view` or `getindex` return another stack with the method 
applied to all layers.
"""
struct GeoStack{T,D,R,M} <: AbstractGeoStack
    data::T
    dims::D
    refdims::R
    metadata::M
end
GeoStack(data::Vararg{<:AbstractGeoArray}; 
         keys=Symbol.(name.(data)), dims=(), refdims=(), metadata=nothing) =
    GeoStack(NamedTuple{keys}(data), dims, refdims, metadata)
GeoStack(s::AbstractGeoStack, keys=keys(s)) = begin
    data = convert.(GeoArray, (s[key] for key in keys)) |> NamedTuple{keys}
    GeoStack(data, dims(s), refdims(s), metadata(s))
end

metadata(s::GeoStack) = s.metadata
refdims(s::GeoStack) = s.refdims
rebuild(s::GeoStack, data, newdims, newref) =
    GeoStack(data, newdims, newref, metadata(s))

Base.parent(s::GeoStack) = s.data
Base.getindex(s::GeoStack, key::Symbol) = getindex(parent(s), key)
Base.getindex(s::GeoStack, I...) = 
    rebuild(s, (a[I...] for a in values(s)) |> NamedTuple{keys(s)})
Base.view(s::GeoStack, I...) = 
    rebuild(s, (view(a, I...) for a in values(s)) |> NamedTuple{keys(s)})

DimensionalData.select(s::GeoStack, I, args...) = 
    rebuild(s, (select(a, I, args...) for a in values(s)) |> NamedTuple{keys(s)})
