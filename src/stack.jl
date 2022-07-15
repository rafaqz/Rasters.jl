# Accept either Symbol or String keys, but allways convert to Symbol
const Key = Union{Symbol,AbstractString}

"""
    AbstractRasterStack

Abstract supertype for objects that hold multiple [`AbstractRaster`](@ref)
that share spatial dimensions.

They are `NamedTuple`-like structures that may either contain `NamedTuple`
of [`AbstractRaster`](@ref), string paths that will load [`AbstractRaster`](@ref),
or a single path that points to as a file itself containing multiple layers, like
NetCDF or HDF5. Use and syntax is similar or identical for all cases.

`AbstractRasterStack` can hold layers that share some or all of their dimensions.
They cannot have the same dimension with t different length or spatial extent as
another layer.

`getindex` on a `AbstractRasterStack` generally returns a memory backed standard
[`Raster`](@ref). `geoarray[:somelayer] |> plot` plots the layers array,
while `geoarray[:somelayer, X(1:100), Band(2)] |> plot` will plot the
subset without loading the whole array.

`getindex` on a `AbstractRasterStack` with a key returns another stack with
getindex applied to all the arrays in the stack.
"""
abstract type AbstractRasterStack{L} <: AbstractDimStack{L} end

missingval(stack::AbstractRasterStack) = stack.missingval
filename(stack::AbstractRasterStack) = filename(parent(stack))
missingval(s::AbstractRasterStack, key::Symbol) = _singlemissingval(missingval(s), key)

isdisk(A::AbstractRasterStack) = isdisk(first(A))

setcrs(x::AbstractRasterStack, crs) = set(x, setcrs(dims(x), crs)...)
setmappedcrs(x::AbstractRasterStack, mappedcrs) = set(x, setmappedcrs(dims(x), mappedcrs)...)

"""
    subset(s::AbstractRasterStack, keys)

Subset a stack to hold only the layers in `keys`, where `keys` is a `Tuple`
or `Array` of `String` or `Symbol`, or a `Tuple` or `Array` of `Int`
"""
subset(s::AbstractRasterStack, keys) = subset(s, Tuple(keys))
function subset(s::AbstractRasterStack, keys::NTuple{<:Any,<:Key})
    RasterStack(map(k -> s[k], Tuple(keys)))
end
function subset(s::AbstractRasterStack, I::NTuple{<:Any,Int})
    subset(s, map(i -> keys(s)[i], I))
end

_singlemissingval(mvs::NamedTuple, key) = mvs[key]
_singlemissingval(mv, key) = mv

# DimensionalData methods ######################################################

# Always read a stack before loading it as a table.
DD.DimTable(stack::AbstractRasterStack) = invoke(DD.DimTable, Tuple{AbstractDimStack}, read(stack))

function DD.layers(s::AbstractRasterStack{<:FileStack{<:Any,Keys}}) where Keys
    NamedTuple{Keys}(map(K -> s[K], Keys))
end
function DD.layers(s::AbstractRasterStack{<:OpenStack{<:Any,Keys}}) where Keys
    NamedTuple{Keys}(map(K -> s[K], Keys))
end

function DD.rebuild(
    s::AbstractRasterStack, data, dims=dims(s), refdims=refdims(s), 
    layerdims=DD.layerdims(s), metadata=metadata(s), layermetadata=DD.layermetadata(s),
    missingval=missingval(s), 
)
    DD.basetypeof(s)(data, dims, refdims, layerdims, metadata, layermetadata, missingval)
end
function DD.rebuild(s::AbstractRasterStack;
    data=parent(s), dims=dims(s), refdims=refdims(s), layerdims=DD.layerdims(s),
    metadata=metadata(s), layermetadata=DD.layermetadata(s),
    missingval=missingval(s),
)
    DD.basetypeof(s)(
        data, dims, refdims, layerdims, metadata, layermetadata, missingval
    )
end

function DD.rebuild_from_arrays(
    s::AbstractRasterStack, das::NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractDimArray}}}; 
    refdims=DD.refdims(s), 
    metadata=DD.metadata(s), 
    data=map(parent, das), 
    dims=DD.combinedims(das...), 
    layerdims=map(DD.basedims, das),
    layermetadata=map(DD.metadata, das),
    missingval=map(missingval, das),
)
    rebuild(s; data, dims, refdims, layerdims, metadata, layermetadata, missingval)
end

# Base methods #################################################################

Base.names(s::AbstractRasterStack) = keys(s)
Base.copy(stack::AbstractRasterStack) = map(copy, stack)

#### Stack getindex ####
# Different to DimensionalData as we construct a Raster
Base.getindex(s::AbstractRasterStack, key::AbstractString) = s[Symbol(key)]
function Base.getindex(s::AbstractRasterStack, key::Symbol)
    data_ = parent(s)[key]
    dims_ = dims(s, DD.layerdims(s, key))
    metadata = DD.layermetadata(s, key)
    Raster(data_, dims_, refdims(s), key, metadata, missingval(s, key))
end
# Key + Index
@propagate_inbounds @inline function Base.getindex(s::AbstractRasterStack, key::Symbol, i1, I...)
    A = s[key][i1, I...]
end


# Concrete AbstrackRasterStack implementation #################################################

"""
    RasterStack <: AbstrackRasterStack

    RasterStack(data...; name, kw...)
    RasterStack(data::Union{Vector,Tuple}; name, kw...)
    RasterStack(data::NamedTuple; kw...))
    RasterStack(s::AbstractRasterStack; kw...)
    RasterStack(s::AbstractRaster; layersfrom=Band, kw...)

Load a file path or a `NamedTuple` of paths as a `RasterStack`, or convert arguments, a 
`Vector` or `NamedTuple` of `Raster` to `RasterStack`.

# Arguments

- `data`: A `NamedTuple` of [`Raster`](@ref), or a `Vector`, `Tuple` or splatted arguments
    of [`Raster`](@ref). The latter options must pass a `name` keyword argument.

# Keywords

- `name`: Used as stack layer names when a `Tuple`, `Vector` or splat of `Raster` is passed in.
- `metadata`: A `Dict` or `DimensionalData.Metadata` object.
- `refdims`: `Tuple` of `Dimension` that the stack was sliced from.
- `layersfrom`: `Dimension` to source stack layers from if the file is not already multi-layered.
    `nothing` is default, so that a single `RasterStack(raster)` is a single layered stack.
    `RasterStack(raster; layersfrom=Band)` will use the bands as layers.
- `lazy`: A `Bool` specifying if to load the stack lazily from disk. `false` by default.

```julia
files = (:temp="temp.tif", :pressure="pressure.tif", :relhum="relhum.tif")
stack = RasterStack(files; mappedcrs=EPSG(4326))
stack[:relhum][Lat(Contains(-37), Lon(Contains(144))
```
"""
struct RasterStack{L<:Union{FileStack,OpenStack,NamedTuple},D<:Tuple,R<:Tuple,LD<:NamedTuple,M,LM,MV} <: AbstractRasterStack{L}
    data::L
    dims::D
    refdims::R
    layerdims::LD
    metadata::M
    layermetadata::LM
    missingval::MV
end
# Multi-file stack from strings
function RasterStack(
    filenames::Union{AbstractArray{<:AbstractString},Tuple{<:AbstractString,Vararg}};
    name=map(filekey, filenames), keys=name, kw...
)
    RasterStack(NamedTuple{Tuple(keys)}(Tuple(filenames)); kw...)
end
function RasterStack(filenames::NamedTuple{K,<:Tuple{<:AbstractString,Vararg}};
    crs=nothing, mappedcrs=nothing, source=nothing, lazy=false, kw...
) where K
    layers = map(keys(filenames), values(filenames)) do key, fn
        source = source isa Nothing ? _sourcetype(fn) : source
        crs = defaultcrs(source, crs)
        mappedcrs = defaultmappedcrs(source, mappedcrs)
        _open(source, fn; key) do ds
            data = if lazy
                FileArray(ds, fn; key)
            else
                _open(Array, source, ds; key)
            end
            dims = DD.dims(ds, crs, mappedcrs)
            md = metadata(ds)
            mv = missingval(ds)
            Raster(data, dims; name=key, metadata=md, missingval=mv)
        end
    end
    RasterStack(NamedTuple{K}(layers); kw...)
end
# Multi Raster stack from Tuple of AbstractArray
function RasterStack(data::Tuple{Vararg{<:AbstractArray}}, dims::Tuple; name=nothing, keys=name, kw...)
    return RasterStack(NamedTuple{cleankeys(keys)}(data), dims; kw...)
end
# Multi Raster stack from NamedTuple of AbstractArray
function RasterStack(data::NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractArray}}}, dims::Tuple; kw...)
    # TODO: make this more sophisticated an match dimension length to axes?
    layers = map(data) do A
        Raster(A, dims[1:ndims(A)])
    end
    return RasterStack(layers; kw...)
end
# Multi Raster stack from AbstractDimArray splat
RasterStack(layers::AbstractDimArray...; kw...) = RasterStack(layers; kw...)
# Multi Raster stack from tuple with `keys` keyword
function RasterStack(layers::Tuple{Vararg{<:AbstractRaster}}; 
    name=map(name, layers), keys=name, kw...
)
    RasterStack(NamedTuple{cleankeys(keys)}(layers); kw...)
end
# Multi Raster stack from NamedTuple
function RasterStack(layers::NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractRaster}}};
    resize=nothing, refdims=(), metadata=NoMetadata(), kw...
)
    # resize if not matching sizes - resize can be `crop`, `resample` or `extend`
    layers = resize isa Nothing ? layers : resize(layers)
    # DD.comparedims(layers...)
    dims=DD.combinedims(layers...)
    data=map(parent, layers)
    layerdims=map(DD.basedims, layers)
    layermetadata=map(DD.metadata, layers)
    missingval=map(Rasters.missingval, layers)
    return RasterStack(
        data, dims, refdims, layerdims, metadata,
        layermetadata, missingval
    )
end
# Single-file stack from a string
function RasterStack(filename::AbstractString;
    dims=nothing, refdims=(), metadata=nothing, crs=nothing, mappedcrs=nothing,
    layerdims=nothing, layermetadata=nothing, missingval=nothing,
    source=_sourcetype(filename), name=nothing, keys=name, layersfrom=nothing,
    resize=nothing, lazy=true, ext=nothing
)
    st = if isfile(filename)
        st = if haslayers(_sourcetype(filename))
            crs = defaultcrs(source, crs)
            mappedcrs = defaultmappedcrs(source, mappedcrs)
            data, field_kw = _open(filename) do ds
                dims = dims isa Nothing ? DD.dims(ds, crs, mappedcrs) : dims
                refdims = refdims == () || refdims isa Nothing ? () : refdims
                layerdims = layerdims isa Nothing ? DD.layerdims(ds) : layerdims
                metadata = metadata isa Nothing ? DD.metadata(ds) : metadata
                layermetadata = layermetadata isa Nothing ? DD.layermetadata(ds) : layermetadata
                missingval = missingval isa Nothing ? Rasters.missingval(ds) : missingval
                data = FileStack{source}(ds, filename; keys)
                data, (; dims, refdims, layerdims, metadata, layermetadata, missingval)
            end
            RasterStack(data; field_kw...)
        else
            # Band dims acts as layers
            RasterStack(Raster(filename; lazy); layersfrom)
        end
        # Maybe split the stack into separate arrays to remove extra dims.
        if !(keys isa Nothing)
            map(identity, st)
        else
            st
        end
    elseif isdir(filename)
        # Load a whole directory
        filenames = readdir(filename)
        length(filenames) > 0 || throw(ArgumentError("No files in directory $filename"))
        # Detect keys from names
        keys = if isnothing(keys)
            all_shared = true
            stripped = lstrip.(x -> x in (" ", "_"), (x -> x[1:end]).(filenames))
            Symbol.(replace.(first.(splitext.(stripped)), Ref(" " => "_")))
        else
            keys
        end
        RasterStack(joinpath.(Ref(filename), filenames); keys)
    end
    return lazy ? st : read(st)
end
function RasterStack(A::Raster; 
    layersfrom=nothing, name=nothing, keys=name, metadata=metadata(A), refdims=refdims(A), kw...
)

    keys = keys isa Union{AbstractString,Symbol,Name} ? (keys,) : keys
    layers = if isnothing(layersfrom)
        keys = if keys isa Nothing
            keys = DD.name(A) in (NoName(), Symbol(""), Name(Symbol(""))) ? ("layer1",) : DD.name(A) 
        else
            keys
        end
        NamedTuple{cleankeys(keys)}((A,))
    else
        keys = keys isa Nothing ? _layerkeysfromdim(A, layersfrom) : keys
        slices = slice(A, layersfrom)
        NamedTuple{cleankeys(keys)}(Tuple(slices))
    end
    RasterStack(layers; refdims=refdims, metadata=metadata, kw...)
end
# Stack from stack, dims args
RasterStack(st::AbstractRasterStack, dims::Tuple; kw...) = RasterStack(st; dims, kw...)
# Stack from table, dims args
function RasterStack(table, dims::Tuple; name=_not_a_dimcol(table, dims), keys=name, kw...)
    # TODO use `name` everywhere, not keys
    if keys isa Symbol
        col = Tables.getcolumn(table, keys)
        layers = NamedTuple{(keys,)}((reshape(col, map(length, dims)),))
    else
        layers = map(keys) do k
            col = Tables.getcolumn(table, k)
            reshape(col, map(length, dims))
        end |> NamedTuple{keys}
    end
    RasterStack(layers, dims; kw...)
end

function DD.modify(f, s::AbstractDimStack{<:FileStack})
    open(s) do o 
        map(a -> modify(f, a), o)
    end
end

function Base.open(f::Function, s::AbstractDimStack{<:FileStack})
    f(rebuild(s; data=OpenStack(parent(s))))
end

function _layerkeysfromdim(A, dim)
    map(index(A, dim)) do x
        if x isa Number
            Symbol(string(DD.dim2key(dim), "_", x))
        else
            Symbol(x)
        end
    end
end

# Rebuild from internals
function RasterStack(
    data::Union{FileStack,OpenStack,NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractArray}}}};
    dims, refdims=(), layerdims, metadata=NoMetadata(), layermetadata, missingval) 
    st = RasterStack(
        data, dims, refdims, layerdims, metadata, layermetadata, missingval
    )
end
# RasterStack from another stack
function RasterStack(s::AbstractDimStack; name=cleankeys(Base.keys(s)), keys=name,
    data=NamedTuple{keys}(s[key] for key in keys),
    dims=dims(s), refdims=refdims(s), layerdims=DD.layerdims(s),
    metadata=metadata(s), layermetadata=DD.layermetadata(s),
    missingval=missingval(s)
)
    st = RasterStack(
        data, DD.dims(s), refdims, layerdims, metadata, layermetadata, missingval
    )

    # TODO This is a bit of a hack, it should use `formatdims`. 
    return set(st, dims...)
end

Base.convert(::Type{RasterStack}, src::AbstractDimStack) = RasterStack(src)

Raster(stack::RasterStack) = cat(values(stack)...; dims=Band([keys(stack)...]))

defaultcrs(T::Type, crs) = crs
defaultcrs(T::Type, ::Nothing) = defaultcrs(T)
defaultcrs(T::Type) = nothing
defaultmappedcrs(T::Type, crs) = crs
defaultmappedcrs(T::Type, ::Nothing) = defaultmappedcrs(T)
defaultmappedcrs(T::Type) = nothing

# Precompile
precompile(RasterStack, (String,))

@deprecate stack(args...; kw...) RasterStack(args...; kw...)
