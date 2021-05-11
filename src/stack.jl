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

"""
    childkwargs(s::AbstractGeoStack)

Returns the keyword arguments that will be passed to the child array constructor.
"""
childkwargs(s::AbstractGeoStack) = s.childkwargs

window(stack::AbstractGeoStack) = stack.window


# Interface methods ############################################################

"""
    modify(f, series::AbstractGeoStack)

Apply function `f` to the data of the child `AbstractGeoArray`s.

`f` must return an idenically sized array.

This method triggers a complete rebuild of all objects,
and disk based objects will be transferred to memory.

This is useful for swapping out array backend for an
entire stack to `CuArray` from CUDA.jl to copy data to a GPU,
and potentially other types like `DAarray` from Distributed.jl.
"""
DD.modify(f, s::AbstractGeoStack) = GeoStack(s; data=_mapdata(A -> modify(f, A), s))

# Base methods #################################################################

# @propagate_inbounds function Base.getindex(s::AbstractGeoStack, I...)
    # rebuild(s; data=NamedTuple{keys(s)}(a[I...] for a in values(s)))
# end
# Dict/Array hybrid with dims
# @propagate_inbounds function Base.getindex(s::AbstractGeoStack, key::Key, I::Vararg{<:Dimension})
    # getindex(s, key, DD.dims2indices(dims(s, key), I)...)
# end

Base.names(s::AbstractGeoStack) = keys(s)
Base.copy(stack::AbstractGeoStack) = rebuild(stack; data=map(copy, getsource(stack)))

"""
    Base.cat(stacks::AbstractGeoStack...; [keys=keys(stacks[1])], dims)

Concatenate all or a subset of layers for all passed in stacks.

# Keywords

- `keys`: `Tuple` of `Symbol` for the stack keys to concatenate.
- `dims`: Dimension of child array to concatenate on.

# Example

Concatenate the :sea_surface_temp and :humidity layers in the time dimension:

```julia
cat(stacks...; keys=(:sea_surface_temp, :humidity), dims=Ti)
```
"""
function Base.cat(stacks::AbstractGeoStack...; keys=keys(stacks[1]), dims)
    vals = Tuple(cat((s[key] for s in stacks)...; dims=dims) for key in keys)
    GeoStack(stacks[1], data=NamedTuple{keys}(vals))
end

_mapdata(f, s::AbstractGeoStack) = NamedTuple{cleankeys(keys(s))}(map(f, values(s)))


function DD.rebuild(s::T;
    data=data(s), dims=dims(s), refdims=refdims(s), layerdims=DD.layerdims(s),
    metadata=metadata(s), layermetadata=DD.layermetadata(s), 
    window=window(s), childkwargs=childkwargs(s)
) where T<:AbstractGeoStack
    DD.basetypeof(T)(data, dims, refdims, layerdims, metadata, layermetadata, window, childkwargs)
end

# getsource(s::MemGeoStack{<:NamedTuple}, args...) = data(s, args...)

childdata(f, childobj, ::AbstractGeoStack) = f(childobj)


#### Stack getindex ####
# Symbol key
@propagate_inbounds function Base.getindex(s::AbstractGeoStack, key::Symbol) 
    GeoArray(
        data(s)[key], dims(s, DD.layerdims(s, key)), refdims(s), key, DD.layermetadata(s, key), missing
    )
end
@propagate_inbounds function Base.getindex(s::AbstractGeoStack, key::Symbol, i1, I...) 
    readwindowed(s[key], window(s), i1, I...)
end
@propagate_inbounds function Base.getindex(s::AbstractDimStack, i1::Int, I::Int...)
    map(A -> Base.getindex(A, i1, I...), data(s))
end

# @propagate_inbounds function Base.view(s::AbstractGeoStack, I...)
#     rebuild(s; data=NamedTuple{keys(s)}(view(a, I...) for a in values(s)))
# end

# @propagate_inbounds function Base.getindex(s::AbstractGeoStack, key::Symbol)
#     A = data(s)[key]
#     window_ = maybewindow2indices(dims(s, key), window(s))
#     readwindowed(A, window_)
#     GeoArray(File
# end
# @propagate_inbounds function Base.getindex(
#     s::AbstractGeoStack, key::Symbol, i1::StandardIndices, I::StandardIndices...
# )
#     A = rebuild(data(s)[key]; childkwargs(s)...)
#     window_ = maybewindow2indices(dims(A, key), window(s))
#     readwindowed(A, window_, i1, I...)
# end


# Concrete AbstrackGeoStack implementation ######################################################

"""
    GeoStack <: AbstrackGeoStack

    GeoStack(data...; keys, kwargs...)
    GeoStack(data::Union{Vector,Tuple}; keys, kwargs...)
    GeoStack(data::NamedTuple; window=(), metadata=NoMetadata(), refdims=(), childkwargs=()) =
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
- `childkwargs`: A `NamedTuple` of keyword arguments to pass to the constructor.
- `refdims`: `Tuple` of  position `Dimension` the array was sliced from.
"""
struct GeoStack{L,D,R,LD,M,LM,W,K} <: AbstractGeoStack{L}
    data::L
    dims::D
    refdims::R
    layerdims::LD
    metadata::M
    layermetadata::LM
    window::W
    childkwargs::K
end
GeoStack(das::AbstractDimArray...; kw...) = GeoStack(das; kw...)
function GeoStack(filenames::Uniont{AbstractArray,Tuple}; keys=map(filekey, filenames), kw...)
    stack(NamedTuple{Tuple(keys)}(Tuple(filenames)); kw...)
end
function GeoStack(data::Tuple{Vararg{<:AbstractGeoArray}}; keys=map(name, data), kw...)
    GeoStack(NamedTuple{cleankeys(keys)}(data); kw...)
end
function GeoStack(das::NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractGeoArray}}}; 
    data=map(parent, das), dims=DD.combinedims(das...), refdims=(), 
    layerdims=map(basedims, das), metadata=NoMetadata(), 
    layermetadata=map(DD.metadata, das), window=(), childkwargs=(),
)
    GeoStack(data, dims, refdims, layerdims, metadata, layermetadata, window, childkwargs)
end
function GeoStack(data::FileStack; 
    dims, refdims=(), layerdims, metadata=NoMetadata(), layermetadata=(), window=(), childkwargs=(),
)
    GeoStack(data, dims, refdims, layerdims, metadata, layermetadata, window, childkwargs)
end
function GeoStack(s::AbstractDimStack; keys=cleankeys(Base.keys(s)),
    data=NamedTuple{keys}(s[key] for key in keys),
    dims=dims(s), refdims=refdims(s), layerdims=DD.layerdims(s), 
    window=(), metadata=metadata(s), layermetadata=DD.layermetadata(s), childkwargs=()
)
    GeoStack(data, dims, refdims, layerdims, metadata, layermetadata, window, childkwargs)
end
function GeoStack(filenames::NamedTuple; crs=nothing, mappedcrs=nothing, kw...)
    layerfields = map(keys(filenames), values(filenames)) do k, fn
        readdata(fn) do ds
            data = FileArray(ds, fn, k)
            md = metadata(ds)
            dims = DD.dims(ds, crs, mappedcrs)
            (; data, dims, keys, md)
        end
    end
    data = map(f-> f.data, layerfields)
    dims = DD.commondims(map(f-> f.dims, layerfields)...)
    layerdims = map(f-> DD.basedims(f.dims), layerfields)
    layermetadata = map(f-> f.md, layerfields)
    GeoStack(data; dims, layerdims, layermetadata, kw...)
end
function GeoStack(filename::AbstractString; 
    metadata=nothing, window=(), crs=nothing, mappedcrs=nothing
)
    dims, layerdims, keys, metadata, sizes = _ncread(filename) do ds
        keys = Tuple(map(Symbol, _nondimkeys(ds)))
        md = metadata isa Nothing ? DD.metadata(ds) : metadata
        dims = DD.dims(ds, crs, mappedcrs)
        ldims = _ncdbasedims(ds)
        sizes = _ncdsizes(ds, keys)
        dims, ldims, keys, md, sizes
    end
    data = FileStack{_NCD,keys}(filename, sizes) 
    GeoStack(data; 
        dims=dims, layerdims=layerdims, metadata=metadata, 
        layermetadata=map(_->NoMetadata(), layerdims), childkwargs=(),
    )
end

Base.convert(::Type{GeoStack}, src::AbstractDimStack) = GeoStack(src)
