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
abstract type AbstractRasterStack{L} <: DD.AbstractDimStack{L} end

missingval(stack::AbstractRasterStack) = getfield(stack, :missingval)
missingval(s::AbstractRasterStack, name::Symbol) = _singlemissingval(missingval(s), name)
filename(stack::AbstractRasterStack{<:NamedTuple}) = map(s -> filename(s), stack)
filename(stack::AbstractRasterStack{<:Union{FileStack,OpenStack}}) = filename(parent(stack))

isdisk(st::AbstractRasterStack) = isdisk(layers(st, 1))

setcrs(x::AbstractRasterStack, crs) = set(x, setcrs(dims(x), crs)...)
setmappedcrs(x::AbstractRasterStack, mappedcrs) = set(x, setmappedcrs(dims(x), mappedcrs)...)

_singlemissingval(mvs::NamedTuple, name) = mvs[name]
_singlemissingval(mv, name) = mv

# DimensionalData methods ######################################################

# Always read a stack before loading it as a table.
DD.DimTable(stack::AbstractRasterStack) = invoke(DD.DimTable, Tuple{DD.AbstractDimStack}, read(stack))

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
    refdims=refdims(s),
    metadata=DD.metadata(s),
    dims=nothing,
    layerdims=map(DD.basedims, das),
    layermetadata=map(DD.metadata, das),
    missingval=map(missingval, das),
)
    data = map(parent, das)
    if dims isa Nothing
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

DD.name(s::AbstractRasterStack) = keys(s)
Base.names(s::AbstractRasterStack) = keys(s)
Base.copy(stack::AbstractRasterStack) = map(copy, stack)

#### Stack getindex ####
# Different to DimensionalData as we construct a Raster
Base.getindex(s::AbstractRasterStack, name::AbstractString) = s[Symbol(name)]
function Base.getindex(s::AbstractRasterStack, name::Symbol)
    data_ = parent(s)[name]
    dims_ = dims(s, DD.layerdims(s, name))
    metadata = DD.layermetadata(s, name)
    mv = missingval(s, name)
    mv = mv isa eltype(data_) ? mv : nothing
    Raster(data_, dims_, refdims(s), name, metadata, mv)
end

# Concrete AbstrackRasterStack implementation #################################################

"""
    RasterStack <: AbstrackRasterStack

    RasterStack(data...; name, kw...)
    RasterStack(data::Union{Vector,Tuple}; name, kw...)
    RasterStack(data::NamedTuple; kw...))
    RasterStack(data::RasterStack; kw...)
    RasterStack(data::Raster; layersfrom=Band, kw...)
    RasterStack(filepath::AbstractString; kw...)

Load a file path or a `NamedTuple` of paths as a `RasterStack`, or convert arguments, a
`Vector` or `NamedTuple` of `Raster`s to `RasterStack`.

# Arguments

- `data`: A `NamedTuple` of [`Raster`](@ref)s or `String`, or a `Vector`, `Tuple` or splatted
    arguments of [`Raster`](@ref). The latter options must pass a `name` keyword argument.
- `filepath`: A file (such as netcdf or tif) to be loaded as a stack, or a directory path
    containing multiple files.

# Keywords

- `name`: Used as stack layer names when a `Tuple`, `Vector` or splat of `Raster` is passed in.
    Has no effect when `NameTuple` is used - the `NamedTuple` keys are the layer names.
$GROUP_KEYWORD 
- `metadata`: A `Dict` or `DimensionalData.Metadata` object.
- `missingval`: a single value for all layers or a `NamedTuple` of
    missingval for each layer. `nothing` specifies no missing value.
$CONSTRUCTOR_CRS_KEYWORD 
$CONSTRUCTOR_MAPPEDCRS_KEYWORD 
- `refdims`: `Tuple` of `Dimension` that the stack was sliced from.

For when one or multiple filepaths are used:

$DROPBAND_KEYWORD
$LAZY_KEYWORD
$SOURCE_KEYWORD

For when a single `Raster` is used:

- `layersfrom`: `Dimension` to source stack layers from if the file is not already multi-layered.
    `nothing` is default, so that a single `RasterStack(raster)` is a single layered stack.
    `RasterStack(raster; layersfrom=Band)` will use the bands as layers.

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
# Rebuild from internals
function RasterStack(
    data::Union{FileStack,OpenStack,NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractArray}}}};
    dims::Tuple,
    refdims::Tuple=(),
    layerdims,
    metadata=NoMetadata(),
    layermetadata,
    missingval,
    crs=nokw,
    mappedcrs=nokw,
    name=nokw,
    resize=nokw,
)
    st = RasterStack(
        data, dims, refdims, layerdims, metadata, layermetadata, missingval
    )
    dims = format(dims, CartesianIndices(st))
    dims = crs isa NoKW ? dims : setcrs(dims, crs)
    dims = mappedcrs isa NoKW ? dims : setmappedcrs(dims, mappedcrs)
    st = rebuild(st; dims)
    return name isa NoKW ? st : st[Dimensions._astuple(name)]
end
# Convert Tuple/Array of array to NamedTuples using name/key
function RasterStack(data::Tuple{Vararg{<:AbstractArray}}, dims::Tuple;
    name::Union{Tuple,AbstractArray,NamedTuple,Nothing}=nothing, 
    kw...
)
    isnothing(name) && throw(ArgumentError("pass a Tuple, Array or NamedTuple of names to the `name` keyword"))
    return RasterStack(NamedTuple{cleankeys(name)}(data); dims, kw...)
end
# Multi Raster stack from NamedTuple of AbstractArray
function RasterStack(data::NamedTuple{<:Any,<:Tuple{Vararg{<:AbstractArray}}}, dims::Tuple; kw...)
    # TODO: make this more sophisticated and match dimension length to axes?
    # We dont worry about Raster keywords because these rasters will be deconstructed 
    # again later, and `kw` will define the RasterStack keywords
    layers = map(data) do A
        Raster(A, dims[1:ndims(A)])
    end
    return RasterStack(layers; kw...)
end
# Multi Raster stack from AbstractDimArray splat
RasterStack(layers::AbstractDimArray...; kw...) = RasterStack(layers; kw...)
# Multi Raster stack from tuple with `name` keyword
function RasterStack(layers::Tuple{Vararg{<:AbstractDimArray}};
    name=map(name, layers),
    kw...
)
    RasterStack(NamedTuple{cleankeys(name)}(layers); kw...)
end
# Multi RasterStack from NamedTuple
# This method is called after most other RasterStack methods.
function RasterStack(layers::NamedTuple{K,<:Tuple{Vararg{<:AbstractDimArray}}};
    resize::Union{Function,NoKW}=nokw,
    _layers=resize isa NoKW ? layers : resize(layers),
    dims::Tuple=DD.combinedims(_layers...),
    refdims::Tuple=(),
    missingval=map(missingval, _layers),
    metadata=NoMetadata(),
    layermetadata=map(DD.metadata, _layers),
    layerdims::NamedTuple{K}=map(DD.basedims, _layers),
    kw...
) where K
    # Handle values that musbe be `NamedTuple`
    layermetadata = if layermetadata isa NamedTuple
        layermetadata
    elseif layermetadata isa Union{Nothing,NoMetadata}
        map(_ -> NoMetadata(), layers)
    else
        throw(ArgumentError("$layermetadata is not a valid input for `layermetadata`. Try a `NamedTuple` of `Dict`, `MetaData` or `NoMetadata`"))
    end
    missingval = missingval isa NoKW ? nothing : missingval
    data = map(parent, _layers)
    return RasterStack(
        data; dims, refdims, layerdims, metadata, layermetadata, missingval, kw...
    )
end
# Stack from table and dims args
function RasterStack(table, dims::Tuple;
    name=_not_a_dimcol(table, dims),
    kw...
)
    Tables.istable(table) || throw(ArgumentError("object $(typeof(table)) is not a valid input to `RasterStack`"))
    if name isa Symbol
        col = Tables.getcolumn(table, name)
        layers = NamedTuple{(name,)}((reshape(col, map(length, dims)),))
    else
        layers = map(name) do k
            col = Tables.getcolumn(table, k)
            reshape(col, map(length, dims))
        end |> NamedTuple{name}
    end
    return RasterStack(layers; dims, kw...)
end
# Stack from a Raster
function RasterStack(A::AbstractDimArray;
    layersfrom=nokw,
    name=nokw,
    refdims::Tuple=refdims(A),
    metadata=metadata(A),
    missingval=missingval(A),
    kw...
)
    name = name isa Union{AbstractString,Symbol,Name} ? (name,) : name
    layers = if layersfrom isa NoKW
        name = if name isa NoKW
            name = DD.name(A) in (NoName(), Symbol(""), Name(Symbol(""))) ? ("layer1",) : DD.name(A)
        else
            name
        end
        NamedTuple{cleankeys(name)}((A,))
    else
        name = name isa NoKW ? _layerkeysfromdim(A, layersfrom) : name
        slices = slice(A, layersfrom)
        NamedTuple{cleankeys(name)}(slices)
    end
    return RasterStack(layers; refdims, metadata, missingval, kw...)
end
# Stack from stack and dims args
RasterStack(st::DD.AbstractDimStack, dims::Tuple; kw...) = RasterStack(st; dims, kw...)
# RasterStack from another stack
function RasterStack(s::DD.AbstractDimStack;
    data::NamedTuple=parent(s),
    dims::Union{Tuple,NoKW}=dims(s),
    refdims::Tuple=refdims(s),
    layerdims=DD.layerdims(s),
    metadata=metadata(s),
    layermetadata=DD.layermetadata(s),
    missingval=missingval(s),
    kw...
)
    return RasterStack(
        data; dims, refdims, layerdims, metadata, layermetadata, missingval, kw...
    )
end
# Multi-file stack from strings
function RasterStack(
    filenames::Union{AbstractArray{<:AbstractString},Tuple{<:AbstractString,Vararg}};
    name=map(filekey, filenames),
    kw...
)
    RasterStack(NamedTuple{cleankeys(Tuple(name))}(filenames); kw...)
end
function RasterStack(filenames::NamedTuple{K,<:Tuple{<:AbstractString,Vararg}};
    lazy=false,
    dropband=true,
    source=nothing,
    crs=nokw,
    mappedcrs=nokw,
    missingval=nokw,
    kw...
) where K
    layers = map(keys(filenames), values(filenames)) do name, fn
        source = _sourcetrait(fn, source)
        _open(source, fn; name) do ds
            dims = _dims(ds, crs, mappedcrs)
            prod(map(length, dims))
            data = if lazy
                FileArray{typeof(source)}(ds, fn; name)
            else
                _open(source, ds; name) do A
                    _checkmem(A)
                    Array(A)
                end
            end
            md = _metadata(ds)
            missingval = missingval isa NoKW ? Rasters.missingval(ds) : missingval
            raster = Raster(data; dims, name, metadata=md, missingval)
            return dropband ? _drop_single_band(raster, lazy) : raster
        end
    end
    # Try to keep the passed-in missingval as a single value
    RasterStack(NamedTuple{K}(layers); missingval, kw...)
end
# Stack from a String
function RasterStack(filename::AbstractString;
    lazy::Bool=false,
    dropband::Bool=true,
    source::Union{Symbol,Source,NoKW}=nokw,
    name=nokw,
    group=nokw,
    kw...
)
    source = _sourcetrait(filename, source)
    st = if isdir(filename)
        # Load as a whole directory
        filenames = readdir(filename)
        length(filenames) > 0 || throw(ArgumentError("No files in directory $filename"))
        # Detect keys from names
        name = if name isa NoKW
            all_shared = true
            stripped = lstrip.(x -> x in (" ", "_"), (x -> x[1:end]).(filenames))
            Symbol.(replace.(first.(splitext.(stripped)), Ref(" " => "_")))
        else
            name
        end
        RasterStack(joinpath.(Ref(filename), filenames); lazy, kw...)
    else
        # Load as a single file
        if haslayers(source)
            # With multiple named layers
            l_st = _layer_stack(filename; source, name, lazy, group, kw...)

            # Maybe split the stack into separate arrays to remove extra dims.
            if !(name isa NoKW)
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
_open_layers(f, st) = _open_layers(f, st, DD.layers(st), NamedTuple())
function _open_layers(f, st, unopened::NamedTuple{K}, opened::NamedTuple) where K
    open(first(unopened)) do open_layer
        layer_nt = NamedTuple{(first(K),)}((open_layer,))
        _open_layers(f, st, Base.tail(unopened), merge(opened, layer_nt))
    end
end
function _open_layers(f, st, unopened::NamedTuple{()}, opened::NamedTuple)
    f(rebuild(st; data=opened))
end

function _layer_stack(filename;
    source=nokw,
    dims=nokw,
    refdims=(),
    name=nokw,
    group=nokw,
    metadata=nokw,
    layerdims=nokw,
    layermetadata=nokw,
    missingval=nokw,
    crs=nokw,
    mappedcrs=nokw,
    lazy=false,
    kw...
)
    data, field_kw = _open(filename; source) do ds
        layers = _layers(ds, name, group)
        # Create a Dict of dimkey => Dimension to use in `dim` and `layerdims`
        dimdict = _dimdict(ds, crs, mappedcrs)
        refdims = refdims == () || refdims isa Nothing ? () : refdims
        metadata = metadata isa NoKW ? _metadata(ds) : metadata
        layerdims = layerdims isa NoKW ? _layerdims(ds; layers, dimdict) : layerdims
        dims = _sort_by_layerdims(dims isa NoKW ? _dims(ds, dimdict) : dims, layerdims)
        layermetadata = layermetadata isa NoKW ? _layermetadata(ds; layers) : layermetadata
        missingval = missingval isa NoKW ? Rasters.missingval(ds) : missingval
        names = Tuple(map(Symbol, layers.names))
        data = if lazy
            FileStack{typeof(source)}(ds, filename; name=names, group, vars=Tuple(layers.vars))
        else
            arrays = map(layers.vars) do v
                A = Array(v)
                # Hack for NCDatasets.jl bug with zero dimensional arrays
                A isa Array ? A : fill(A)
            end
            NamedTuple{names}(arrays)
        end
        data, (; dims, refdims, layerdims=NamedTuple{names}(layerdims), metadata, layermetadata=NamedTuple{names}(layermetadata), missingval)
    end
    return RasterStack(data; field_kw..., kw...)
end

# Try to sort the dimensions by layer dimension into a sensible
# order that applies without permutation, preferencing the layers
# with most dimensions, and those that come first.
# Intentionally not type-stable
function _sort_by_layerdims(dims, layerdims)
    dimlist = union(layerdims)
    currentorder = nothing
    for i in length(dims):-1:1
        for ldims in dimlist
            length(ldims) == i || continue
            currentorder = _merge_dimorder(ldims, currentorder)
        end
    end
    return DD.dims(dims, currentorder)
end

_merge_dimorder(neworder, ::Nothing) = neworder
function _merge_dimorder(neworder, currentorder)
    ods = otherdims(neworder, currentorder)
    outorder = currentorder
    for od in ods
        # Get the dims position in current order
        i = findfirst(d -> d == od, neworder)
        found = false
        # Find the next dimension that is in the outorder
        for j in 1:length(ods)
            if length(neworder) >= (i + j)
                nextd = neworder[i + j]
                if nextd in outorder
                    n = dimnum(outorder, nextd)
                    outorder = (outorder[1:n-1]..., od, outorder[n:end]...)
                    found = true
                    break
                end
            end
        end
        if !found
            outorder = (outorder..., od)
        end
    end
    return outorder
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
defaultcrs(s::Source, ::NoKW) = defaultcrs(s)
defaultcrs(x) = nothing
defaultmappedcrs(::Source, crs) = crs
defaultmappedcrs(s::Source, ::NoKW) = defaultmappedcrs(s)
defaultmappedcrs(::Source) = nothing
