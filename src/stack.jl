# Accept either Symbol or String keys, but allways convert to Symbol
const Key = Union{Symbol,AbstractString}

const MAX_STACK_SIZE = 200

"""
    AbstractRasterStack

Abstract supertype for objects that hold multiple [`AbstractRaster`](@ref)s
that share spatial dimensions.

They are `NamedTuple`-like structures that may either contain `NamedTuple`
of [`AbstractRaster`](@ref)s, string paths that will load [`AbstractRaster`](@ref)s,
or a single path that points to a file containing multiple layers, like
NetCDF or HDF5. Use and syntax is similar or identical for all cases.

`AbstractRasterStack` can hold layers that share some or all of their dimensions.
They cannot have the same dimension with different length or spatial extent as
another layer.

`getindex` on an `AbstractRasterStack` generally returns a memory backed standard
[`Raster`](@ref). `raster[:somelayer] |> plot` plots the layers array,
while `raster[:somelayer, X(1:100), Band(2)] |> plot` will plot the
subset without loading the whole array.

`getindex` on an `AbstractRasterStack` with a key returns another stack with
`getindex` applied to all the arrays in the stack.
"""
abstract type AbstractRasterStack{L} <: AbstractDimStack{L} end

missingval(stack::AbstractRasterStack) = getfield(stack, :missingval)
filename(stack::AbstractRasterStack) = filename(parent(stack))
missingval(s::AbstractRasterStack, key::Symbol) = _singlemissingval(missingval(s), key)

isdisk(st::AbstractRasterStack) = isdisk(layers(st, 1))

setcrs(x::AbstractRasterStack, crs) = set(x, setcrs(dims(x), crs)...)
setmappedcrs(x::AbstractRasterStack, mappedcrs) = set(x, setmappedcrs(dims(x), mappedcrs)...)

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
    s::AbstractRasterStack{<:Union{FileStack{<:Any,Keys},OpenStack{<:Any,Keys}}}, das::Tuple{Vararg{<:AbstractDimArray}}; kw...
) where Keys
    DD.rebuild_from_arrays(s, NamedTuple{Keys}(das); kw...)
end
function DD.rebuild_from_arrays(
    s::AbstractRasterStack, das::NamedTuple{<:Any,<:Tuple{Vararg{AbstractDimArray}}};
    data=map(parent, das),
    refdims=refdims(s),
    metadata=DD.metadata(s),
    dims=nothing,
    layerdims=map(DD.basedims, das),
    layermetadata=map(DD.metadata, das),
    missingval=map(missingval, das),
)
    if isnothing(dims)
        # invokelatest avoids compiling this for other paths
        Base.invokelatest() do
            dims = DD.combinedims(collect(das))
        end
        rebuild(s; data, dims, refdims, layerdims, metadata, layermetadata, missingval)
    else
        rebuild(s; data, dims, refdims, layerdims, metadata, layermetadata, missingval)
    end
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
    RasterStack(filename::AbstractString; kw...)

Load a file path or a `NamedTuple` of paths as a `RasterStack`, or convert arguments, a
`Vector` or `NamedTuple` of `Raster`s to `RasterStack`.

# Arguments

- `data`: A `NamedTuple` of [`Raster`](@ref)s, or a `Vector`, `Tuple` or splatted arguments
    of [`Raster`](@ref). The latter options must pass a `name` keyword argument.
- `filename`: A file (such as netcdf or tif) to be loaded as a stack, or a directory path
    containing multiple files.

# Keywords

- `name`: Used as stack layer names when a `Tuple`, `Vector` or splat of `Raster` is passed in.
    Has no effect when `NameTuple` is used - the `NamedTuple` keys are the layer names.
- `metadata`: A `Dict` or `DimensionalData.Metadata` object.
- `refdims`: `Tuple` of `Dimension` that the stack was sliced from.
- `layersfrom`: `Dimension` to source stack layers from if the file is not already multi-layered.
    `nothing` is default, so that a single `RasterStack(raster)` is a single layered stack.
    `RasterStack(raster; layersfrom=Band)` will use the bands as layers.
- `lazy`: A `Bool` specifying whether to load the stack lazily from disk. `false` by default.
- `dropband`: drop single band dimensions when creating stacks from filenames. `true` by default.

```julia
files = (temp="temp.tif", pressure="pressure.tif", relhum="relhum.tif")
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
    RasterStack(NamedTuple{cleankeys(Tuple(keys))}(filenames); kw...)
end
function RasterStack(filenames::NamedTuple{K,<:Tuple{<:AbstractString,Vararg}};
    crs=nothing, mappedcrs=nothing, source=nothing, lazy=false, dropband=true, kw...
) where K
    layers = map(keys(filenames), values(filenames)) do key, fn
        source = _sourcetrait(fn, source)
        crs = defaultcrs(source, crs)
        mappedcrs = defaultmappedcrs(source, mappedcrs)
        _open(source, fn; key) do ds
            dims = _dims(ds, crs, mappedcrs)
            prod(map(length, dims))
            data = if lazy
                FileArray{typeof(source)}(ds, fn; key)
            else
                _open(source, ds; key) do A
                    _checkmem(A)
                    Array(A)
                end
            end
            md = _metadata(ds)
            mv = missingval(ds)
            raster = Raster(data, dims; name=key, metadata=md, missingval=mv) 
            return dropband ? _drop_single_band(raster, lazy) : raster
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
    # TODO: make this more sophisticated and match dimension length to axes?
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
        data, dims, refdims, layerdims, metadata, layermetadata, missingval
    )
end
# Stack from a String
function RasterStack(filename::AbstractString;
    source=nothing, name=nothing, keys=name, lazy=false, dropband=true, kw...
)
    source = _sourcetrait(filename, source)
    st = if isdir(filename)
        # Load as a whole directory
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
        RasterStack(joinpath.(Ref(filename), filenames); lazy, kw...)
    else
        # Load as a single file
        if haslayers(source)
            # With multiple named layers
            l_st = _layer_stack(filename; source, name, keys, lazy, kw...)

            # Maybe split the stack into separate arrays to remove extra dims.
            if !(keys isa Nothing)
                map(identity, l_st)
            else
                l_st
            end
        else
            # With bands actings as layers
            RasterStack(Raster(filename; source, lazy, dropband=false); kw...)
        end
    end

    # Maybe drop the Band dimension
    if dropband && hasdim(st, Band()) && size(st, Band()) == 1
         if lazy
             return view(st, Band(1)) # TODO fix dropdims in DiskArrays
         else
             return dropdims(st; dims=Band())
         end
    else
         return st
    end
end

function _layer_stack(filename;
    dims=nothing, refdims=(), metadata=nothing, crs=nothing, mappedcrs=nothing,
    layerdims=nothing, layermetadata=nothing, missingval=nothing,
    source=nothing, name=nothing, keys=name, resize=nothing, lazy=false, kw...
)
    crs = defaultcrs(source, crs)
    mappedcrs = defaultmappedcrs(source, mappedcrs)
    data, field_kw = _open(filename; source) do ds
        layers = _layers(ds, keys)
        # Create a Dict of dimkey => Dimension to use in `dim` and `layerdims`
        dimdict = _dimdict(ds, crs, mappedcrs)
        dims = dims isa Nothing ? _dims(ds, dimdict) : dims
        refdims = refdims == () || refdims isa Nothing ? () : refdims
        metadata = metadata isa Nothing ? _metadata(ds) : metadata
        layerdims = layerdims isa Nothing ? _layerdims(ds; layers, dimdict) : layerdims
        layermetadata = layermetadata isa Nothing ? _layermetadata(ds; layers) : layermetadata
        missingval = missingval isa Nothing ? Rasters.missingval(ds) : missingval
        tuplekeys = Tuple(map(Symbol, layers.keys))
        data = if lazy
            FileStack{typeof(source)}(ds, filename; keys=tuplekeys, vars=Tuple(layers.vars))
        else
            NamedTuple{tuplekeys}(map(Array, layers.vars))
        end
        data, (; dims, refdims, layerdims=NamedTuple{tuplekeys}(layerdims), metadata, layermetadata=NamedTuple{tuplekeys}(layermetadata), missingval)
    end
    return RasterStack(data; field_kw..., kw...)
end

# Stack from a Raster
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
    return RasterStack(layers; refdims=refdims, metadata=metadata, kw...)
end
# Stack from stack and dims args
RasterStack(st::AbstractRasterStack, dims::Tuple; kw...) = RasterStack(st; dims, kw...)
# Stack from table and dims args
function RasterStack(table, dims::Tuple; name=_not_a_dimcol(table, dims), keys=name, kw...)
    Tables.istable(table) || throw(ArgumentError("object $(typeof(table)) is not a valid input to `RasterStack`"))
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
    return RasterStack(layers, dims; kw...)
end
# Rebuild from internals
function RasterStack(
    data::Union{FileStack,OpenStack,NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractArray}}}};
    dims, refdims=(), layerdims, metadata=NoMetadata(), layermetadata, missingval)
    return RasterStack(
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

function DD.modify(f, s::AbstractRasterStack{<:FileStack{<:Any,K}}) where K
    open(s) do o
        map(K) do k
            Array(parent(ost)[k])
        end
    end
end

# Open a single file stack
function Base.open(f::Function, st::AbstractRasterStack{<:FileStack{<:Any,K}}; kw...) where K
    ost = OpenStack(parent(st))
    layers = map(K) do k
        ost[k]
    end |> NamedTuple{K}
    out = f(rebuild(st; data=layers))
    close(ost)
    return out
end
# Open a multi-file stack or just apply f to a memory backed stack
function Base.open(f::Function, st::AbstractRasterStack{<:NamedTuple}; kw...)
    isdisk(st) ? _open_layers(f, st) : f(st)
end

# Open all layers through nested closures, applying `f` to the rebuilt open stack
_open_layers(f, st) = _open_layers(f, st, DD.layers(f), NamedTuple())
function _open_layers(f, st, unopened::NamedTuple{K}, opened::NamedTuple) where K
    open(first(unopened)) do open_layer
        _open_layers(f, st, Base.tail(unopened), merge(opened, NamedTuple{(first(K))}(open_layer)))
    end
end
function _open_layers(f, st, unopened::NamedTuple{()}, opened)
    f(rebuild(st; data=opened))
end

function _layerkeysfromdim(A, dim)
    hasdim(A, dim) || throw(ArgumentError("`layersrom` dim `$(dim2key(dim))` not found in `$(map(basetypeof, dims(A)))`"))
    vals = parent(lookup(A, dim))
    l = length(vals)
    if l > MAX_STACK_SIZE 
        D = basetypeof(dim)
        options = map(basetypeof, dims(otherdims(A, dim), d -> length(d) <= MAX_STACK_SIZE))
        toolong = "Lookup of `layersfrom=$D` is too long to use for stack layers: $l"
        if length(options) > 0
            throw(ArgumentError("$toolong. Choose a different dimension from: $options"))
        else
            throw(ArgumentError("$toolong. Maybe try a simple `Raster` ?"))
        end
    end
    map(vals) do x
        if x isa Number
            Symbol(string(DD.dim2key(dim), "_", x))
        else
            Symbol(x)
        end
    end
end

Base.convert(::Type{RasterStack}, src::AbstractDimStack) = RasterStack(src)

Raster(stack::RasterStack) = cat(values(stack)...; dims=Band([keys(stack)...]))

defaultcrs(::Source, crs) = crs
defaultcrs(s::Source, ::Nothing) = defaultcrs(s)
defaultcrs(x) = nothing
defaultmappedcrs(::Source, crs) = crs
defaultmappedcrs(s::Source, ::Nothing) = defaultmappedcrs(s)
defaultmappedcrs(::Source) = nothing
