# Accept either Symbol or String keys, but allways convert to Symbol
const Key = Union{Symbol,AbstractString}

"""
    AbstractGeoStack

Abstract supertype for objects that hold multiple [`AbstractGeoArray`](@ref)
that share spatial bounds.

They are `NamedTuple`-like structures that may either contain `NamedTuple`
of [`AbstractGeoArray`](@ref), string paths that will load [`AbstractGeoArray`](@ref),
or a single path that points to as a file itself containing multiple layers, like
NetCDF or HDF5. Use and syntax is similar or identical for all cases.

`getindex` on a `AbstractGeoStack` generally returns a memory backed standard
[`GeoArray`](@ref). `geoarray[:somelayer] |> plot` plots the layers array,
while `geoarray[:somelayer, X(1:100), Band(2)] |> plot` will plot the
subset without loading the whole array.

`getindex` on a `AbstractGeoStack` with a key returns another stack with
getindex applied to all the arrays in the stack.
"""
abstract type AbstractGeoStack{L} <: AbstractDimStack{L} end

layermissingval(stack::AbstractGeoStack) = stack.layermissingval
filename(stack::AbstractGeoStack) = filename(data(stack))
missingval(s::AbstractGeoStack, key::Symbol) = _singlemissingval(layermissingval(s), key)

_singlemissingval(mvs::NamedTuple, key) = mvs[key]
_singlemissingval(mv, key) = mv

# DimensionalData methods ######################################################

# Always read a stack before loading it as a table.
DD.DimTable(stack::AbstractGeoStack) = invoke(DD.DimTable, Tuple{AbstractDimStack}, read(stack))

function DD.layers(s::AbstractGeoStack{<:FileStack{<:Any,Keys}}) where Keys
    NamedTuple{Keys}(map(K -> s[K], Keys))
end

function DD.rebuild(s::T;
    data=data(s), dims=dims(s), refdims=refdims(s), layerdims=DD.layerdims(s),
    metadata=metadata(s), layermetadata=DD.layermetadata(s),
    layermissingval=layermissingval(s), 
) where T<:AbstractGeoStack
    DD.basetypeof(T)(
        data, dims, refdims, layerdims, metadata, layermetadata, layermissingval
    )
end

function DD.rebuild_from_arrays(
    s::AbstractGeoStack, das::NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractDimArray}}}; 
    refdims=DD.refdims(s), 
    metadata=DD.metadata(s), 
    data=map(parent, das), 
    dims=DD.combinedims(das...), 
    layerdims=map(DD.basedims, das),
    layermetadata=map(DD.metadata, das),
    layermissingval=map(missingval, das),
)
    rebuild(s; data, dims, refdims, layerdims, metadata, layermetadata, layermissingval)
end

# Base methods #################################################################

Base.names(s::AbstractGeoStack) = keys(s)
Base.copy(stack::AbstractGeoStack) = map(copy, stack)

#### Stack getindex ####
# Different to DimensionalData as we construct a GeoArray
@propagate_inbounds function Base.getindex(s::AbstractGeoStack, key::Symbol)
    data_ = data(s)[key]
    dims_ = dims(s, DD.layerdims(s, key))
    metadata = DD.layermetadata(s, key)
    GeoArray(data_, dims_, refdims(s), key, metadata, missingval(s, key))
end
# Key + Index
@propagate_inbounds @inline function Base.getindex(s::AbstractGeoStack, key::Symbol, i1, I...)
    A = s[key][i1, I...]
end


# Concrete AbstrackGeoStack implementation #################################################

"""
    GeoStack <: AbstrackGeoStack

    GeoStack(data...; keys, kwargs...)
    GeoStack(data::Union{Vector,Tuple}; keys, kwargs...)
    GeoStack(data::NamedTuple; window=nothing, metadata=NoMetadata(), refdims=()))
    GeoStack(s::AbstractGeoStack; [keys, data, refdims, window, metadata])

A concrete `AbstractGeoStack` implementation. Holds layers of [`GeoArray`](@ref).

# Arguments

- `data`: A `NamedTuple` of [`GeoArray`](@ref), or a `Vector`, `Tuple` or splatted arguments
    of [`GeoArray`](@ref). The latter options must pass a `keys` keyword argument.

# Keywords

- `keys`: Used as stack keys when a `Tuple` or `Vector` or splat of geoarrays are passed in.
- `window`: A `Tuple` of `Dimension`/`Selector`/indices that will be applied to the
    contained arrays when they are accessed.
- `refdims`: Reference dimensions from earlier subsetting.
- `metadata`: A `DimensionalData.Metadata` object.
- `refdims`: `Tuple` of  position `Dimension` the array was sliced from.
"""
struct GeoStack{L,D,R,LD,M,LM,LMV} <: AbstractGeoStack{L}
    data::L
    dims::D
    refdims::R
    layerdims::LD
    metadata::M
    layermetadata::LM
    layermissingval::LMV
end
# Multi-file stack from strings
function GeoStack(
    filenames::Union{AbstractArray{<:AbstractString},Tuple{<:AbstractString,Vararg}};
    keys=map(filekey, filenames), kw...
)
    GeoStack(NamedTuple{Tuple(keys)}(Tuple(filenames)); kw...)
end
function GeoStack(filenames::NamedTuple{K,<:Tuple{<:AbstractString,Vararg}};
    crs=nothing, mappedcrs=nothing, source=nothing, kw...
) where K
    layers = map(keys(filenames), values(filenames)) do key, fn
        source = source isa Nothing ? _sourcetype(fn) : source
        crs = defaultcrs(source, crs)
        mappecrs = defaultmappedcrs(source, mappedcrs)
        _open(fn; key) do ds
            data = FileArray(ds, fn; key)
            dims = DD.dims(ds, crs, mappedcrs)
            md = metadata(ds)
            mv = missingval(ds)
            GeoArray(data, dims; name=key, metadata=md, missingval=mv)
        end
    end
    GeoStack(NamedTuple{K}(layers); kw...)
end
# Multi GeoArray stack from splat with `keys` keyword
GeoStack(layers::AbstractDimArray...; kw...) = GeoStack(layers; kw...)
# Multi GeoArray stack from tuple with `keys` keyword
function GeoStack(layers::Tuple{Vararg{<:AbstractGeoArray}}; keys=map(name, layers), kw...)
    GeoStack(NamedTuple{cleankeys(keys)}(layers); kw...)
end
# Multi GeoArray stack from NamedTuple
function GeoStack(layers::NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractGeoArray}}};
    resize=nothing, refdims=(), metadata=NoMetadata(), window=nothing, kw...
)
    # resize if not matching sizes - resize can be `crop` or `extent`
    layers = resize isa Nothing ? layers : resize(layers)
    # DD.comparedims(layers...)
    dims=DD.combinedims(layers...)
    data=map(parent, layers)
    layerdims=map(DD.basedims, layers)
    layermetadata=map(DD.metadata, layers)
    layermissingval=map(missingval, layers)
    st = GeoStack(
        data, dims, refdims, layerdims, metadata,
        layermetadata, layermissingval
    )
    # For NamedTuple stacks we apply the whole window
    window === nothing ? st : view(st, window...)
end
# Single-file stack from a string
function GeoStack(filename::AbstractString;
    dims=nothing, refdims=(), metadata=nothing, crs=nothing, mappedcrs=nothing,
    layerdims=nothing, layermetadata=nothing, layermissingval=nothing,
    source=_sourcetype(filename), keys=nothing, window=nothing
)
    crs = defaultcrs(source, crs)
    mappedcrs = defaultmappedcrs(source, mappedcrs)
    data, field_kw = _open(filename) do ds
        dims = dims isa Nothing ? DD.dims(ds, crs, mappedcrs) : dims
        refdims = refdims == () || refdims isa Nothing ? DD.refdims(ds, filename) : refdims
        layerdims = layerdims isa Nothing ? DD.layerdims(ds) : layerdims
        metadata = metadata isa Nothing ? DD.metadata(ds) : metadata
        layermetadata = layermetadata isa Nothing ? DD.layermetadata(ds) : layermetadata
        layermissingval = layermissingval isa Nothing ? GeoData.layermissingval(ds) : layermissingval
        data = FileStack{source}(ds, filename; keys)
        data, (; dims, refdims, layerdims, metadata, layermetadata, layermissingval)
    end
    GeoStack(data; field_kw..., window)
end



# Rebuild from internals
function GeoStack(
    data::Union{FileStack,NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractArray}}}};
    dims, refdims=(), layerdims, metadata=NoMetadata(), layermetadata, layermissingval, window=nothing) 
    st = GeoStack(
        data, dims, refdims, layerdims, metadata, layermetadata, layermissingval
    )
    window === nothing ? st : view(st, window...)
end
# GeoStack from another stack
function GeoStack(s::AbstractDimStack; keys=cleankeys(Base.keys(s)),
    data=NamedTuple{keys}(s[key] for key in keys),
    dims=dims(s), refdims=refdims(s), layerdims=DD.layerdims(s),
    metadata=metadata(s), layermetadata=DD.layermetadata(s),
    layermissingval=layermissingval(s), window=nothing
)
    st = GeoStack(
        data, dims, refdims, layerdims, metadata, layermetadata, layermissingval
    )
    window === nothing ? st : view(st, window...)
end

Base.convert(::Type{GeoStack}, src::AbstractDimStack) = GeoStack(src)

defaultcrs(T::Type, crs) = crs
defaultcrs(T::Type, ::Nothing) = defaultcrs(T)
defaultcrs(T::Type) = nothing
defaultmappedcrs(T::Type, crs) = crs
defaultmappedcrs(T::Type, ::Nothing) = defaultmappedcrs(T)
defaultmappedcrs(T::Type) = nothing
