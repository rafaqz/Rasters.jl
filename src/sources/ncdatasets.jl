import .NCDatasets
const NCD = NCDatasets

export NCDarray, NCDstack, NCDdimMetadata, NCDarrayMetadata, NCDstackMetadata

"""
    NCDdimMetadata(val::Union{Dict,NamedTuple})
    NCDdimMetadata(pairs::Pair...) => NCDdimMetadata{Dict}
    NCDdimMetadata(; kw...) => NCDdimMetadata{NamedTuple}

`Metadata` wrapper for [`NCDarray`](@ref) dimensions.
"""
struct NCDdimMetadata{T} <: AbstractDimMetadata{T}
    val::T
end

"""
    NCDarrayMetadata(val::Union{Dict,NamedTuple})
    NCDarrayMetadata(pairs::Pair...) => NCDarrayMetadata{Dict}
    NCDarrayMetadata(; kw...) => NCDarrayMetadata{NamedTuple}

`Metadata` wrapper for [`NCDarray`](@ref) metadata.
"""
struct NCDarrayMetadata{T} <: AbstractArrayMetadata{T}
    val::T
end

"""
    NCDstackMetadata(val::Union{Dict,NamedTuple})
    NCDstackMetadata(pairs::Pair...) => NCDstackMetadata{Dict}
    NCDstackMetadata(; kw...) => NCDstackMetadata{NamedTuple}

`Metadata` wrapper for [`NCDarray`](@ref) metadata.
"""
struct NCDstackMetadata{T} <: AbstractStackMetadata{T}
    val::T
end

const UNNAMED_NCD_KEY = "unnamed"

const NCD_FILL_TYPES = (Int8,UInt8,Int16,UInt16,Int32,UInt32,Int64,UInt64,Float32,Float64,Char,String)

# CF standards don't enforce dimension names.
# But these are common, and should take care of most dims.
const DIMMAP = Dict("lat" => Lat,
                    "latitude" => Lat,
                    "lon" => Lon,
                    "long" => Lon,
                    "longitude" => Lon,
                    "time" => Ti,
                    "lev" => Vert,
                    "mlev" => Vert,
                    "level" => Vert,
                    "vertical" => Vert,
                    "x" => X,
                    "y" => Y,
                    "z" => Z,
                   )

# Array ########################################################################
"""
    NCDarray(filename::AbstractString; name=nothing, refdims=(),
             dims=nothing, metadata=nothing, crs=nothing, mappedcrs=EPSG(4326))

A [`DiskGeoArray`](@ref) that loads that loads NetCDF files lazily from disk.

The first non-dimension layer of the file will be used as the array. Dims are usually
detected as [`Lat`](@ref), [`Lon`](@ref), [`Ti`]($DDtidocs), and [`Vert`] or
possibly `X`, `Y`, `Z` when detected. Undetected dims will use the generic `Dim{:name}`.

This is an incomplete implementation of the NetCDF standard. It will currently
handle simple files in lattitude/longitude projections, or projected formats
if you manually specify `crs` and `mappedcrs`. How this is done may also change in
future, including detecting and converting the native NetCDF projection format.

## Arguments

- `filename`: `String` pointing to a netcdf file.

## Keyword arguments

- `crs`: defaults to lat/lon `EPSG(4326)` but may be any `GeoFormat` like `WellKnownText`.
  If the underlying data is in a different projection `crs` will need to be set to allow
  `write` to a different file format. In future this may be detected automatically.
- `dimcrs`: The crs projection actually present in the Dimension index `Vector`, which
  may be different to the underlying projection. Defaults to lat/lon `EPSG(4326)` but
  may be any crs `GeoFormat`.
- `name`: `Symbol` name for the array. Will use array key if not supplied.
- `dims`: `Tuple` of `Dimension`s for the array. Detected automatically, but can be passed in.
- `refdims`: `Tuple of` position `Dimension`s the array was sliced from.
- `missingval`: Value reprsenting missing values. Detected automatically when possible, but
  can be passed it.
- `metadata`: [`Metadata`](@ref) object for the array. Detected automatically as
  [`NCDarrayMetadata`](@ref), but can be passed in.

## Example

```julia
A = NCDarray("folder/file.ncd")
# Select Australia from the default lat/lon coords:
A[Lat(Between(-10, -43), Lon(Between(113, 153)))
```
"""
struct NCDarray{T,N,A,D<:Tuple,R<:Tuple,Na<:Symbol,Me,Mi,S,K
               } <: DiskGeoArray{T,N,D,LazyArray{T,N}}
    filename::A
    dims::D
    refdims::R
    name::Na
    metadata::Me
    missingval::Mi
    size::S
    key::K
end
function NCDarray(filename::AbstractString, key...; kw...)
    isfile(filename) || error("File not found: $filename")
    _ncread(dataset -> NCDarray(dataset, filename, key...; kw...), filename)
end
# Safe file-loading wrapper method. We always open the datset and close
# it again when we are done.
function NCDarray(dataset::NCD.Dataset, filename, key=nothing;
    crs=nothing,
    mappedcrs=EPSG(4326),
    name=nothing,
    dims=nothing,
    refdims=(),
    metadata=nothing,
    missingval=missing,
)
    keys_ = _nondimkeys(dataset)
    key = (key isa Nothing) ? first(keys_) : key |> Symbol
    var = dataset[string(key)]
    dims = dims isa Nothing ? GeoData.dims(dataset, key, crs, mappedcrs) : dims
    name = Symbol(name isa Nothing ? key : name)
    metadata_ = metadata isa Nothing ? GeoData.metadata(var, GeoData.metadata(dataset)) : metadata
    size_ = map(length, dims)
    T = eltype(var)
    N = length(dims)

    NCDarray{T,N,typeof.((filename,dims,refdims,name,metadata_,missingval,size_,key))...
       }(filename, dims, refdims, name, metadata_, missingval, size_, key)
end

key(A::NCDarray) = A.key

# AbstractGeoArray methods

withsourcedata(f, A::NCDarray) = _ncread(ds -> f(ds[string(key(A))]), filename(A))

# Base methods

Base.size(A::NCDarray) = A.size

"""
    Base.write(filename::AbstractString, ::Type{NCDarray}, s::AbstractGeoArray)

Write an NCDarray to a NetCDF file using NCDatasets.jl

Returns `filename`.
"""
function Base.write(filename::AbstractString, ::Type{NCDarray}, A::AbstractGeoArray)
    dataset = NCD.Dataset(filename, "c")
    try
        println("    Writing netcdf...")
        _ncwritevar!(dataset, A)
    finally
        close(dataset)
    end
    return filename
end

# Stack ########################################################################


"""
    NCDstack(filename::String; refdims=(), window=(), metadata=nothing, childkwargs=())
    NCDstack(filenames; keys, kw...)
    NCDstack(filenames...; keys, kw...)
    NCDstack(filenames::NamedTuple; refdims=(), window=(), metadata=nothing, childkwargs=())

A lazy [`AbstractGeoStack`](@ref) that uses NCDatasets.jl to load NetCDF files.
Can load a single multi-layer netcdf file, or multiple single-layer netcdf
files. In multi-file mode it returns a regular `GeoStack` with a `childtype`
of [`NCDarray`](@ref).

Indexing into `NCDstack` with layer keys (`Symbol`s) returns a [`GeoArray`](@ref).
Dimensions are usually detected as [`Lat`](@ref), [`Lon`](@ref), [`Ti`]($DDtidocs),
and [`Vert`] or `X`, `Y`, `Z` when detected. Undetected dims use the generic `Dim{:name}`.

# Arguments

- `filename`: `Tuple` or `Vector` or splatted arguments of `String`,
  or single `String` path, to NetCDF files.

# Keyword arguments

- `keys`: Used as stack keys when a `Tuple`, `Vector` or splat of filenames are passed in.
  These default to the first non-dimension data key in each NetCDF file.
- `refdims`: Add dimension position array was sliced from. Mostly used programatically.
- `window`: A `Tuple` of `Dimension`/`Selector`/indices that will be applied to the
  contained arrays when they are accessed.
- `metadata`: A [`StackMetadata`](@ref) object.
- `childkwargs`: A `NamedTuple` of keyword arguments to pass to the
  [`NCDarray`](@ref) constructor.

# Examples

```julia
stack = NCDstack(filename; window=(Lat(Between(20, 40),))
# Or
stack = NCDstack([fn1, fn1, fn3, fn4])
# And index with a layer key
stack[:soiltemp]
```
"""
struct NCDstack{T,R,W,M,K} <: DiskGeoStack{T}
    filename::T
    refdims::R
    window::W
    metadata::M
    childkwargs::K
end
function NCDstack(filename::AbstractString;
    refdims=(),
    window=(),
    metadata=_ncread(metadata, filename),
    childkwargs=()
)
    NCDstack(filename, refdims, window, metadata, childkwargs)
end
# These actually return a DiskStack
NCDstack(filenames::AbstractString...; kw...) = NCDstack(filenames; kw...)
function NCDstack(filenames; keys=_ncfilenamekeys(filenames), kw...)
    DiskStack(filenames; keys=keys, childtype=NCDarray, kw...)
end
NCDstack(filenames::NamedTuple; kw...) = DiskStack(filenames; childtype=NCDarray, kw...)

childtype(::NCDstack) = NCDarray
childkwargs(stack::NCDstack) = stack.childkwargs
crs(stack::NCDarray) = get(childkwargs(stack), :crs, EPSG(4326))
mappedcrs(stack::NCDarray) = get(childkwargs(stack), :mappedcrs, nothing)
missingval(::NCDstack) = missing

# AbstractGeoStack methods
withsource(f, ::Type{NCDarray}, path::AbstractString, key=nothing) = _ncread(f, path)
withsourcedata(f, ::Type{NCDarray}, path::AbstractString, key) =
    _ncread(d -> f(d[string(key)]), path)

# Override the default to get the dims of the specific key,
# and pass the crs and dimscrs from the stack
DD.dims(s::NCDstack, dataset, key::Key) = dims(dataset, key, crs(s), mappedcrs(s))

# Base methods

Base.keys(s::NCDstack{<:AbstractString}) = cleankeys(_ncread(_nondimkeys, getsource(s)))

"""
    Base.write(filename::AbstractString, ::Type{NCDstack}, s::AbstractGeoStack)

Write an NCDstack to a single netcdf file, using NCDatasets.jl.

Currently [`DimMetadata`](@ref) is not handled, and [`ArrayMetadata`](@ref)
from other [`AbstractGeoArray`](@ref) @types is ignored.
"""
function Base.write(filename::AbstractString, ::Type{<:NCDstack}, s::AbstractGeoStack)
    dataset = NCD.Dataset(filename, "c"; attrib=_stringdict(metadata(s)))
    try
        map(key -> _ncwritevar!(dataset, s[key]), keys(s))
    finally
        close(dataset)
    end
end

# DimensionalData methods for NCDatasets types ###############################

function DD.dims(dataset::NCD.Dataset, key::Key, crs=nothing, mappedcrs=nothing)
    v = dataset[string(key)]
    dims = []
    for (i, dimname) in enumerate(NCD.dimnames(v))
        if haskey(dataset, dimname)
            dvar = dataset[dimname]
            # Find the matching dimension constructor. If its an unknown name use
            # the generic Dim with the dim name as type parameter
            dimtype = haskey(DIMMAP, dimname) ? DIMMAP[dimname] : DD.basetypeof(DD.key2dim(Symbol(dimname)))
            index = dvar[:]
            meta = NCDdimMetadata(DD.metadatadict(dvar.attrib))
            mode = _ncdmode(index, dimtype, crs, mappedcrs, meta)

            # Add the dim containing the dimension var array
            push!(dims, dimtype(index, mode, meta))
        else
            # The var doesn't exist. Maybe its `complex` or some other marker,
            # so make it a custom `Dim` with `NoIndex`
            push!(dims, Dim{Symbol(dimname)}(1:size(v, i), NoIndex(), nothing))
        end
    end
    (dims...,)
end

DD.metadata(dataset::NCD.Dataset) = NCDstackMetadata(DD.metadatadict(dataset.attrib))
DD.metadata(dataset::NCD.Dataset, key::Key) = metadata(dataset[string(key)])
DD.metadata(var::NCD.CFVariable) = NCDarrayMetadata(DD.metadatadict(var.attrib))
DD.metadata(var::NCD.CFVariable, stackmetadata::NCDstackMetadata) = begin
    md = NCDarrayMetadata(DD.metadatadict(var.attrib))
    md[:_stack] = stackmetadata
    md
end

missingval(var::NCD.CFVariable) = missing

# Direct loading: better memory handling?
# readwindowed(A::NCD.CFVariable) = readwindowed(A, axes(A)...)
# readwindowed(A::NCD.CFVariable, i, I...) = begin
#     var = A.var
#     indices = to_indices(var, (i, I...))
#     shape = Base.index_shape(indices...)
#     dest = Array{eltype(var),length(shape)}(undef, map(length, shape)...)
#     NCD.load!(var, dest, indices...)
#     dest
# end

# Utils ########################################################################

function _ncfilenamekeys(filenames)
    cleankeys(_ncread(ds -> first(_nondimkeys(ds)), fn) for fn in filenames)
end

function _ncdmode(index::AbstractArray{<:Number}, dimtype, crs, mappedcrs, metadata)
    # Assume the locus is at the center of the cell if boundaries aren't provided.
    # http://cfconventions.org/cf-conventions/cf-conventions.html#cell-boundaries
    # Unless its a time dimension.
    order = _ncdorder(index)
    span = _ncdspan(index, order)
    sampling = Intervals(Center())
    if dimtype in (Lat, Lon)
        # If the index is regularly spaced and there is no crs
        # then there is probably just one crs - the mappedcrs
        crs = if crs isa Nothing && span isa Regular
            mappedcrs
        else
            crs
        end
        Mapped(order, span, sampling, crs, mappedcrs)
    else
        Sampled(order, span, sampling)
    end
end
function _ncdmode(index::AbstractArray{<:Dates.AbstractTime}, dimtype, crs, mappedcrs, metadata)
    order = _ncdorder(index)
    span, sampling  = _get_period(index, metadata)
    Sampled(order, span, sampling)
end
_ncdmode(index, dimtype, crs, mappedcrs, mode) = Categorical()

function _ncdorder(index)
    index[end] > index[1] ? Ordered(ForwardIndex(), ForwardArray(), ForwardRelation()) :
                            Ordered(ReverseIndex(), ReverseArray(), ForwardRelation())
end

function _ncdspan(index, order)
    # Handle a length 1 index
    length(index) == 1 && return Regular(zero(eltype(index)))
    step = index[2] - index[1]
    for i in 2:length(index) -1
        # If any step sizes don't match, its Irregular
        if !(index[i+1] - index[i] ≈ step)
            bounds = if length(index) > 1
                beginhalfcell = abs((index[2] - index[1]) * 0.5)
                endhalfcell = abs((index[end] - index[end-1]) * 0.5)
                if DD.isrev(indexorder(order))
                    index[end] - endhalfcell, index[1] + beginhalfcell
                else
                    index[1] - beginhalfcell, index[end] + endhalfcell
                end
            else
                index[1], index[1]
            end
            return Irregular(bounds)
        end
    end
    # Otherwise regular
    return Regular(step)
end

function _get_period(index, metadata::NCDdimMetadata)
    if haskey(metadata, :delta_t)
        period = _parse_period(metadata[:delta_t])
        period isa Nothing || return Regular(period), Points()
    elseif haskey(metadata, :avg_period)
        period = _nc_parse_period(metadata[:avg_period])
        period isa Nothing || return Regular(period), Intervals(Center())
    end
    return sampling = Irregular(), Points()
end

function _parse_period(period_str::String)
    regex = r"(\d\d\d\d)-(\d\d)-(\d\d) (\d\d):(\d\d):(\d\d)"
    mtch = match(regex, period_str)
    if mtch isa Nothing
        nothing
    else
        vals = Tuple(parse.(Int, mtch.captures))
        periods = (Year, Month, Day, Hour, Minute, Second)
        if length(vals) == length(periods)
            sum(map((p, v) -> p(v), periods, vals))
        end
    end
end

_stringdict(metadata) = attrib = Dict(string(k) => v for (k, v) in metadata)

_ncread(f, path::String) = NCD.Dataset(f, path)

function _nondimkeys(dataset)
    dimkeys = keys(dataset.dim)
    removekeys = if "bnds" in dimkeys
        dimkeys = setdiff(dimkeys, ("bnds",))
        boundskeys = [dataset[k].attrib["bounds"] for k in dimkeys if haskey(dataset[k].attrib, "bounds")]
        union(dimkeys, boundskeys)
    else
        dimkeys
    end
    setdiff(keys(dataset), removekeys)
end

# Add a var array to a dataset before writing it.
function _ncwritevar!(dataset, A::AbstractGeoArray{T,N}) where {T,N}
    A = reorder(A, ForwardIndex()) |> a -> reorder(a, ForwardRelation())
    # Define required dim vars
    for dim in dims(A)
        key = lowercase(string(name(dim)))
        haskey(dataset.dim, key) && continue

        # Shift index before conversion to Mapped
        dim = _ncshiftindex(dim)
        if dim isa Lat || dim isa Lon
            dim = convertmode(Mapped, dim)
        end

        md = metadata(dim)
        # TODO handle dim attribs
        attribvec = [] #md isa Nothing ? [] : [val(md)...]
        NCD.defDim(dataset, key, length(dim))
        println("        key: \"", key, "\" of type: ", eltype(dim))
        NCD.defVar(dataset, key, Vector(index(dim)), (key,); attrib=attribvec)
    end
    # TODO actually convert the metadata types
    attrib = if metadata isa NCDarrayMetadata
        _stringdict(metadata(A))
    else
        Dict()
    end
    # Remove stack metdata if it is attached
    pop!(attrib, :_stack, nothing)
    # Set _FillValue
    if ismissing(missingval(A))
        eltyp = _notmissingtype(Base.uniontypes(T)...)
        fillval = NCD.fillvalue(eltyp)
        A = replace_missing(A, fillval)
        attrib["_FillValue"] = fillval
    elseif missingval(A) isa T
        attrib["_FillValue"] = missingval(A)
    else
        @warn "`missingval` $(missingval(A)) is not the same type as your data $T."
    end

    key = if string(name(A)) == ""
        UNNAMED_NCD_KEY
    else
        string(name(A))
    end
    println("        key: \"", key, "\" of type: ", T)
    dimnames = lowercase.(string.(name.(dims(A))))
    attribvec = [attrib...]
    var = NCD.defVar(dataset, key, eltype(A), dimnames; attrib=attribvec)

    var[:] = data(A)
end

_notmissingtype(::Type{Missing}, next...) = _notmissingtype(next...)
_notmissingtype(x::Type, next...) = x in NCD_FILL_TYPES ? x : _notmissingtype(next...)
_notmissingtype() = error("Your data is not a type that netcdf can store")

_ncshiftindex(dim::Dimension) = _ncshiftindex(mode(dim), dim)
_ncshiftindex(::IndexMode, dim::Dimension) = dim
function _ncshiftindex(mode::AbstractSampled, dim::Dimension)
    if span(mode) isa Regular && sampling(mode) isa Intervals
        # We cant easily shift a DateTime value
        if eltype(dim) isa Dates.AbstractDateTime
            if !(loci(dim) isa Center)
                @warn "To save to netcdf, DateTime values should be the interval Center, rather than the $(nameof(typeof(loci(dim))))"
            end
            dim
        else
            shiftindexloci(Center(), dim)
        end
    else
        dim
    end
end
