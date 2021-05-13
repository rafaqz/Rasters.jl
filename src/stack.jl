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

window(stack::AbstractGeoStack) = stack.window
layermissingval(stack::AbstractGeoStack) = stack.layermissingval
missingval(s::AbstractGeoStack, key::Symbol) = _singlemissingval(layermissingval(s), key)
filename(stack::AbstractGeoStack) = filename(data(stack))

_singlemissingval(mvs::NamedTuple, key) = mvs[key]
_singlemissingval(mv, key) = mv

# Base methods #################################################################

Base.names(s::AbstractGeoStack) = keys(s)
Base.copy(stack::AbstractGeoStack) = rebuild(stack; data=map(copy, stack))

function DD.layers(s::AbstractGeoStack{<:FileStack{<:Any,Keys}}) where Keys
    NamedTuple{Keys}(map(K -> s[K], Keys))
end

function DD.rebuild(s::T;
    data=data(s), dims=dims(s), refdims=refdims(s), layerdims=DD.layerdims(s),
    metadata=metadata(s), layermetadata=DD.layermetadata(s), 
    layermissingval=layermissingval(s), window=window(s),
) where T<:AbstractGeoStack
    DD.basetypeof(T)(
        data, dims, refdims, layerdims, metadata, layermetadata, layermissingval, window
    )
end

#### Stack getindex ####
# Symbol key
@propagate_inbounds function Base.getindex(s::AbstractGeoStack, key::Symbol) 
    data_ = data(s)[key]
    dims_ = dims(s, DD.layerdims(s, key))
    A = GeoArray(data_, dims_, refdims(s), key, DD.layermetadata(s, key), missingval(s, key))
    window(s) == () ? A : view(A, window(s)...)
end
@propagate_inbounds function Base.getindex(s::AbstractGeoStack, key::Symbol, i1, I...) 
    s[key][i1, I...]
end
@propagate_inbounds function Base.getindex(s::AbstractGeoStack, i1::Int, I::Int...)
    window(s) == () ? A[i1, I...] : map(A -> view(A, window(s)...)[i1, I...], data(s))
end

# @propagate_inbounds function Base.view(s::AbstractGeoStack, I...)
#     rebuild(s; data=NamedTuple{keys(s)}(view(a, I...) for a in values(s)))
# end

# Concrete AbstrackGeoStack implementation ######################################################

"""
    GeoStack <: AbstrackGeoStack

    GeoStack(data...; keys, kwargs...)
    GeoStack(data::Union{Vector,Tuple}; keys, kwargs...)
    GeoStack(data::NamedTuple; window=(), metadata=NoMetadata(), refdims=()))
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
struct GeoStack{L,D,R,LD,M,LM,LMV,W} <: AbstractGeoStack{L}
    data::L
    dims::D
    refdims::R
    layerdims::LD
    metadata::M
    layermetadata::LM
    layermissingval::LMV
    window::W
end
GeoStack(das::AbstractDimArray...; kw...) = GeoStack(das; kw...)
function GeoStack(data::Tuple{Vararg{<:AbstractGeoArray}}; keys=map(name, data), kw...)
    GeoStack(NamedTuple{cleankeys(keys)}(data); kw...)
end
function GeoStack(das::NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractGeoArray}}}; 
    data=map(parent, das), dims=DD.combinedims(das...), refdims=(), 
    layerdims=map(DD.basedims, das), metadata=NoMetadata(), 
    layermetadata=map(DD.metadata, das), 
    layermissingval=map(missingval, das), window=(),
)
    GeoStack(
        data, dims, refdims, layerdims, metadata, layermetadata, 
        layermissingval, window,
    )
end
function GeoStack(data::NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractArray}}}; 
    dims, refdims=(), layerdims, metadata=NoMetadata(), layermetadata, 
    layermissingval, window=(),
)
    GeoStack(
        data, dims, refdims, layerdims, metadata, layermetadata, 
        layermissingval, window,
    )
end
function GeoStack(data::FileStack; 
    dims, refdims=(), layerdims, metadata=NoMetadata(), layermetadata=(), 
    layermissingval, window=(),
)
    GeoStack(
        data, dims, refdims, layerdims, metadata, layermetadata, 
        layermissingval, window,
    )
end
function GeoStack(s::AbstractDimStack; keys=cleankeys(Base.keys(s)),
    data=NamedTuple{keys}(s[key] for key in keys),
    dims=dims(s), refdims=refdims(s), layerdims=DD.layerdims(s), 
    metadata=metadata(s), layermetadata=DD.layermetadata(s), 
    layermissingval=layermissingval(s), window=()
)
    GeoStack(data, dims, refdims, layerdims, metadata, layermetadata, layermissingval, window)
end
function GeoStack(
    filenames::Union{AbstractArray{<:AbstractString},Tuple{<:AbstractString,Vararg}}; 
    keys=map(filekey, filenames), kw...
)
    @show keys filenames kw
    GeoStack(NamedTuple{Tuple(keys)}(Tuple(filenames)); kw...)
end
# Multi-file stack from strings
function GeoStack(filenames::NamedTuple{K,<:Tuple{<:AbstractString,Vararg}}; 
    crs=nothing, mappedcrs=nothing, kw...
) where K
    layerfields = map(keys(filenames), values(filenames)) do k, fn
        source = _sourcetype(fn)
        crs = defaultcrs(source, crs)
        mappecrs = defaultmappedcrs(source, mappedcrs)
        _read(fn, k) do ds
            data = FileArray(ds, fn, k)
            md = metadata(ds)
            dims = DD.dims(ds, crs, mappedcrs)
            mv = missingval(ds)
            (; data, dims, keys, md, mv)
        end
    end
    layerfields = NamedTuple{K}(layerfields)
    data = map(f-> f.data, layerfields)
    @show map(f-> f.dims, layerfields)
    # TODO this should be uniondims(
    dims = DD.commondims(map(f-> f.dims, layerfields)...)
    layerdims = map(f-> DD.basedims(f.dims), layerfields)
    layermetadata = map(f-> f.md, layerfields)
    layermissingval = map(f-> f.mv, layerfields)
    GeoStack(data; dims, layerdims, layermetadata, layermissingval, kw...)
end
# Single-file stack from a string
function GeoStack(filename::AbstractString; 
    refdims=(), metadata=nothing, crs=nothing, mappedcrs=nothing, kw...
)
    source = _sourcetype(filename)
    crs = defaultcrs(source, crs)
    mappedcrs = defaultmappedcrs(source, mappedcrs)
    data, field_kw = _read(filename) do ds
        keys = Tuple(map(Symbol, layerkeys(ds)))
        dims = DD.dims(ds, crs, mappedcrs)
        refdims = refdims == () ? DD.refdims(ds, filename) : refdims
        layerdims = DD.layerdims(ds)
        metadata = metadata isa Nothing ? DD.metadata(ds) : metadata
        layermetadata = DD.layermetadata(ds)
        layermissingval = GeoData.layermissingval(ds)
        sizes = layersizes(ds, keys)
        data = FileStack{source,keys}(filename, sizes) 
        data, (; dims, refdims, layerdims, metadata, layermetadata, layermissingval)
    end
    GeoStack(data; field_kw..., kw...)
end

defaultcrs(T::Type, crs) = crs
defaultmappedcrs(T::Type, crs) = crs
defaultcrs(T::Type, ::Nothing) = defaultcrs(T)
defaultmappedcrs(T::Type, ::Nothing) = defaultmappedcrs(T)
defaultcrs(T::Type) = nothing
defaultmappedcrs(T::Type) = nothing

Base.convert(::Type{GeoStack}, src::AbstractDimStack) = GeoStack(src)

DD.maybestack(As::NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractGeoArray}}}) = GeoStack(As)
