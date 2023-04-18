const CDM = CommonDataModel
const NCD = NCDatasets
const GDS = GRIBDatasets

const CDM_DIM_MAP = Dict(
    "lat" => Y,
    "latitude" => Y,
    "lon" => X,
    "long" => X,
    "longitude" => X,
    "time" => Ti,
    "lev" => Z,
    "mlev" => Z,
    "level" => Z,
    "vertical" => Z,
    "x" => X,
    "y" => Y,
    "z" => Z,
    "band" => Band,
)

const CDM_AXIS_MAP = Dict(
    "X" => X,
    "Y" => Y,
    "Z" => Z,
    "T" => Ti,
)

const CDM_STANDARD_NAME_MAP = Dict(
    "longitude" => X,
    "latitude" => Y,
    "depth" => Z,
    "time" => Ti,
)

_dataset(var::AbstractVariable) = CDM.dataset(var)

haslayers(::Type{<:CDMsource}) = true
defaultcrs(::Type{<:CDMsource}) = EPSG(4326)
defaultmappedcrs(::Type{<:CDMsource}) = EPSG(4326)

# Raster ########################################################################

function Raster(ds::AbstractDataset, filename::AbstractString, key=nothing; kw...)
    if isnothing(key)
        # Find the first valid variable
        for key in layerkeys(ds)
            if ndims(ds[key]) > 0
                @info "No `name` or `key` keyword provided, using first valid layer with name `:$key`"
                return Raster(ds[key], filename, key; source=NCDsource, kw...)
            end
        end
        throw(ArgumentError("dataset at $filename has no array variables"))
    else
       return Raster(ds[key], filename, key; kw...)
    end
end

_firstkey(ds::AbstractDataset, key::Nothing=nothing) = Symbol(first(layerkeys(ds)))
_firstkey(ds::AbstractDataset, key) = Symbol(key)

function FileArray(var::AbstractVariable, filename::AbstractString; kw...)
    da = RasterDiskArray{CDMsource}(var)
    size_ = size(da)
    eachchunk = DA.eachchunk(da)
    haschunks = DA.haschunks(da)
    T = eltype(var)
    N = length(size_)
    FileArray{CDMsource,T,N}(filename, size_; eachchunk, haschunks, kw...)
end

function Base.open(f::Function, A::FileArray{CDMsource}; write=A.write, kw...)
    _open(CDMsource, filename(A); key=key(A), write, kw...) do var
        f(RasterDiskArray{CDMsource}(var, DA.eachchunk(A), DA.haschunks(A)))
    end
end

"""
    Base.write(filename::AbstractString, ::Type{<:CDMsource}, A::AbstractRaster)

Write an NCDarray to a NetCDF file using NCDatasets.jl

Returns `filename`.
"""
function Base.write(filename::AbstractString, ::Type{<:CDMsource}, A::AbstractRaster; append=false, kw...)
    mode  = !isfile(filename) || !append ? "c" : "a";
    ds = NCD.Dataset(filename, mode; attrib=_attribdict(metadata(A)))
    try
        _ncdwritevar!(ds, A; kw...)
    finally
        close(ds)
    end
    return filename
end

# Stack ########################################################################

"""
    Base.write(filename::AbstractString, ::Type{NCDsource}, s::AbstractRasterStack; kw...)

Write an NCDstack to a single netcdf file, using NCDatasets.jl.

Currently `Metadata` is not handled for dimensions, and `Metadata` from other
[`AbstractRaster`](@ref) @types is ignored.

# Keywords

Keywords are passed to `NCDatasets.defVar`.

- `append`: If true, the variable of the current Raster will be appended to
    `filename`. Note that the variable of the current Raster should be not exist
    before. If not, you need to set `append = false`. Rasters.jl can not
    overwrite a previous existing variable.
- `fillvalue`: A value filled in the NetCDF file to indicate missing data. It
    will be stored in the `_FillValue` attribute.
- `chunksizes`: Vector integers setting the chunk size. The total size of a
    chunk must be less than 4 GiB.
- `deflatelevel`: Compression level: 0 (default) means no compression and 9
    means maximum compression. Each chunk will be compressed individually.
- `shuffle`: If true, the shuffle filter is activated which can improve the
    compression ratio.
- `checksum`: The checksum method can be `:fletcher32` or `:nochecksum`
    (checksumming is disabled, which is the default)
 - `typename` (string): The name of the NetCDF type required for vlen arrays
    (https://web.archive.org/save/https://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf-c/nc_005fdef_005fvlen.html)
"""
function Base.write(filename::AbstractString, ::Type{<:CDMsource}, s::AbstractRasterStack; append = false, kw...)
    mode  = !isfile(filename) || !append ? "c" : "a";
    ds = NCD.Dataset(filename, mode; attrib=_attribdict(metadata(s)))
    try
        map(key -> _ncdwritevar!(ds, s[key]), keys(s); kw...)
    finally
        close(ds)
    end
    return filename
end

function create(filename, ::Type{<:CDMsource}, T::Union{Type,Tuple}, dims::DimTuple;
    name=:layer1, keys=(name,), layerdims=map(_->dims, keys), missingval=nothing,
    metadata=NoMetadata(), lazy=true, 
)
    types = T isa Tuple ? T : Ref(T)
    missingval = T isa Tuple ? missingval : Ref(missingval)
    # Create layers of zero arrays
    layers = map(layerdims, keys, types, missingval) do lds, key, t, mv
        A = FillArrays.Zeros{t}(map(length, lds))
        Raster(A, dims=lds; name=key, missingval=mv)
    end
    write(filename, CDMsource, Raster(first(layers)))
    return Raster(filename; source=CDMsource, lazy)
end

# DimensionalData methods for NCDatasets types ###############################

function DD.dims(ds::AbstractDataset, crs=nothing, mappedcrs=nothing)
    map(_dimkeys(ds)) do key
        _ncddim(ds, key, crs, mappedcrs)
    end |> Tuple
end
function DD.dims(var::AbstractVariable, crs=nothing, mappedcrs=nothing)
    names = CDM.dimnames(var)
    map(names) do name
        _ncddim(_dataset(var), name, crs, mappedcrs)
    end |> Tuple
end

_attrib(ds::Union{AbstractDataset, AbstractVariable}) = CDM.attribs(ds)
DD.metadata(ds::AbstractDataset) = _metadatadict(CDMsource, _attrib(ds))
DD.metadata(var::CFVariable) = _metadatadict(CDMsource, _attrib(var))
DD.metadata(var::AbstractVariable) = _metadatadict(CDMsource, _attrib(var))

function DD.layerdims(ds::AbstractDataset)
    keys = Tuple(layerkeys(ds))
    dimtypes = map(keys) do key
        DD.layerdims(ds[string(key)])
    end
    NamedTuple{map(Symbol, keys)}(dimtypes)
end
function DD.layerdims(var::AbstractVariable)
    ds = _dataset(var)
    map(CDM.dimnames(var)) do dimname
        _ncddimtype(_attrib(ds[dimname]), dimname)()
    end |> Tuple
end

DD.layermetadata(ds::AbstractDataset) = _layermetadata(ds, Tuple(layerkeys(ds)))
function _layermetadata(ds, keys)
    dimtypes = map(k -> DD.metadata(ds[string(k)]), keys)
    NamedTuple{map(Symbol, keys)}(dimtypes)
end

missingval(var::CDM.CFVariable{T}) where T = missing isa T ? missing : nothing

function layerkeys(ds::AbstractDataset)
    dimkeys = _dimkeys(ds)
    toremove = if "bnds" in dimkeys
        dimkeys = setdiff(dimkeys, ("bnds",))
        boundskeys = String[]
        for k in dimkeys
            var = ds[k]
            if haskey(_attrib(var), "bounds")
                push!(boundskeys, CDM.attrib(var, "bounds"))
            end
        end
        union(dimkeys, boundskeys)::Vector{String}
    else
        dimkeys::Vector{String}
    end
    return setdiff(keys(ds), toremove)
end

function FileStack(source::Type{<:CDMsource}, ds::AbstractDataset, filename::AbstractString; write, keys)
    keys = map(Symbol, keys isa Nothing ? layerkeys(ds) : keys) |> Tuple
    type_size_ec_hc = map(keys) do key
        var = ds[string(key)]
        Union{Missing,eltype(var)}, size(var), _ncd_eachchunk(var), _ncd_haschunks(var)
    end
    layertypes = map(x->x[1], type_size_ec_hc)
    layersizes = map(x->x[2], type_size_ec_hc)
    eachchunk = map(x->x[3], type_size_ec_hc)
    haschunks = map(x->x[4], type_size_ec_hc)
    return FileStack{source,keys}(filename, layertypes, layersizes, eachchunk, haschunks, write)
end

FileStack{NCDsource}(ds::AbstractDataset, filename::AbstractString; write=false, keys) = FileStack(NCDsource, ds, filename; write, keys)
FileStack{GRIBsource}(ds::AbstractDataset, filename::AbstractString; write=false, keys) = FileStack(GRIBsource, ds, filename; write, keys)

function OpenStack(fs::FileStack{NCDsource,K}) where K
    OpenStack{NCDsource,K}(NCD.Dataset(filename(fs)))
end
function OpenStack(fs::FileStack{GRIBsource,K}) where K
    OpenStack{GRIBsource,K}(GDS.GRIBDataset(filename(fs)))
end
Base.close(os::OpenStack{NCDsource}) = NCD.close(dataset(os))

# In GRIBDatasets, the file is open for reading the values and closed afterwards. 
Base.close(os::OpenStack{GRIBsource}) = nothing

function _open(f, ::Type{NCDsource}, filename::AbstractString; write=false, kw...)
    isfile(filename) || _filenotfound_error(filename)
    mode = write ? "a" : "r"
    NCD.Dataset(filename, mode) do ds
        _open(f, NCDsource, ds; kw...)
    end
end
function _open(f, ::Type{GRIBsource}, filename::AbstractString; write=false, kw...)
    isfile(filename) || _filenotfound_error(filename)
    ds = GRIBDatasets.GRIBDataset(filename)
    _open(f, GRIBsource, ds; kw...)
end

function _open(f, ::Type{<:CDMsource}, ds::AbstractDataset; key=nothing, kw...)
    x = key isa Nothing ? ds : ds[_firstkey(ds, key)]
    cleanreturn(f(x))
end
_open(f, ::Type{<:CDMsource}, var::CDM.CFVariable; kw...) = cleanreturn(f(var))

# Utils ########################################################################

cleanreturn(A::CDM.CFVariable) = Array(A)

# Utils ########################################################################

function _ncddim(ds, dimname::Key, crs=nothing, mappedcrs=nothing)
    if haskey(ds, dimname)
        var = ds[dimname]
        D = _ncddimtype(_attrib(var), dimname)
        lookup = _ncdlookup(ds, dimname, D, crs, mappedcrs)
        return D(lookup)
    else
        # The var doesn't exist. Maybe its `complex` or some other marker,
        # so make it a custom `Dim` with `NoLookup`
        len = _ncfinddimlen(ds, dimname)
        len === nothing && _unuseddimerror()
        lookup = NoLookup(Base.OneTo(len))
        return Dim{Symbol(dimname)}(lookup)
    end
end

function _ncfinddimlen(ds, dimname)
    for key in keys(ds)
        var = ds[key]
        dimnames = CDM.dimnames(var)
        if dimname in dimnames
            return size(var)[findfirst(==(dimname), dimnames)]
        end
    end
    return nothing
end

# Find the matching dimension constructor. If its an unknown name
# use the generic Dim with the dim name as type parameter
function _ncddimtype(attrib, dimname)
    if haskey(attrib, "axis") 
        k = attrib["axis"] 
        if haskey(CDM_AXIS_MAP, k) 
            return CDM_AXIS_MAP[k] 
        end
    end
    if haskey(attrib, "standard_name")
        k = attrib["standard_name"]
        if haskey(CDM_STANDARD_NAME_MAP, k) 
            return CDM_STANDARD_NAME_MAP[k]
        end
    end
    if haskey(CDM_DIM_MAP, dimname) 
        return CDM_DIM_MAP[dimname] 
    end
    return DD.basetypeof(DD.key2dim(Symbol(dimname)))
end

# _ncdlookup
# Generate a `LookupArray` from a netcdf dim.
function _ncdlookup(ds::AbstractDataset, dimname, D::Type, crs, mappedcrs)
    dvar = ds[dimname]
    index = dvar[:]
    metadata = _metadatadict(CDMsource, _attrib(dvar))
    return _ncdlookup(ds, dimname, D, index, metadata, crs, mappedcrs)
end
# For unknown types we just make a Categorical lookup
function _ncdlookup(ds::AbstractDataset, dimname, D::Type, index::AbstractArray, metadata, crs, mappedcrs)
    Categorical(index; order=Unordered(), metadata=metadata)
end
# For Number and AbstractTime we generate order/span/sampling
# We need to include `Missing` in unions in case `_FillValue` is used
# on coordinate variables in a file and propagates here.
function _ncdlookup(
    ds::AbstractDataset, dimname, D::Type, index::AbstractArray{<:Union{Missing,Number,Dates.AbstractTime}},
    metadata, crs, mappedcrs
)
    # Assume the locus is at the center of the cell if boundaries aren't provided.
    # http://cfconventions.org/cf-conventions/cf-conventions.html#cell-boundaries
    order = LA.orderof(index)
    var = ds[dimname]
    if haskey(_attrib(var), "bounds")
        boundskey = var.attrib["bounds"]
        boundsmatrix = Array(ds[boundskey])
        span, sampling = Explicit(boundsmatrix), Intervals(Center())
        return _ncdlookup(D, index, order, span, sampling, metadata, crs, mappedcrs)
    elseif eltype(index) <: Union{Missing,Dates.AbstractTime}
        span, sampling = _ncdperiod(index, metadata)
        return _ncdlookup(D, index, order, span, sampling, metadata, crs, mappedcrs)
    else
        span, sampling = _ncdspan(index, order), Points()
        return _ncdlookup(D, index, order, span, sampling, metadata, crs, mappedcrs)
    end
end
# For X and Y use a Mapped <: AbstractSampled lookup
function _ncdlookup(
    D::Type{<:Union{<:XDim,<:YDim}}, index, order::Order, span, sampling, metadata, crs, mappedcrs
)
    # If the index is regularly spaced and there is no crs
    # then there is probably just one crs - the mappedcrs
    crs = if crs isa Nothing && span isa Regular
        mappedcrs
    else
        crs
    end
    dim = DD.basetypeof(D)()
    return Mapped(index, order, span, sampling, metadata, crs, mappedcrs, dim)
end
# Band dims have a Categorical lookup, with order
function _ncdlookup(D::Type{<:Band}, index, order::Order, span, sampling, metadata, crs, mappedcrs)
    Categorical(index, order, metadata)
end
# Otherwise use a regular Sampled lookup
function _ncdlookup(D::Type, index, order::Order, span, sampling, metadata, crs, mappedcrs)
    Sampled(index, order, span, sampling, metadata)
end

function _ncdspan(index, order)
    # Handle a length 1 index
    length(index) == 1 && return Regular(zero(eltype(index)))
    step = index[2] - index[1]
    for i in 2:length(index)-1
        # If any step sizes don't match, its Irregular
        if !(index[i+1] - index[i] ≈ step)
            bounds = if length(index) > 1
                beginhalfcell = abs((index[2] - index[1]) * 0.5)
                endhalfcell = abs((index[end] - index[end-1]) * 0.5)
                if LA.isrev(order)
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

# delta_t and ave_period are not CF standards, but CDC
function _ncdperiod(index, metadata::Metadata{<:CDMsource})
    if haskey(metadata, "delta_t")
        period = _parse_period(metadata["delta_t"])
        period isa Nothing || return Regular(period), Points()
    elseif haskey(metadata, "avg_period")
        period = _parse_period(metadata["avg_period"])
        period isa Nothing || return Regular(period), Intervals(Center())
    end
    return sampling = Irregular((nothing, nothing)), Points()
end

function _parse_period(period_str::String)
    regex = r"(\d\d\d\d)-(\d\d)-(\d\d) (\d\d):(\d\d):(\d\d)"
    mtch = match(regex, period_str)
    if mtch === nothing
        return nothing
    else
        vals = Tuple(parse.(Int, mtch.captures))
        periods = (Year, Month, Day, Hour, Minute, Second)
        if length(vals) == length(periods)
            compound = sum(map((p, v) -> p(v), periods, vals))
            if length(compound.periods) == 1
                return compound.periods[1]
            else
                return compound
            end
        else
            return nothing
        end
    end
end

_attribdict(md::Metadata{<:CDMsource}) = Dict{String,Any}(string(k) => v for (k, v) in md)
_attribdict(md) = Dict{String,Any}()

_dimkeys(ds::AbstractDataset) = CDM.dimnames(ds)

# Add a var array to a dataset before writing it.
function _ncdwritevar!(ds::AbstractDataset, A::AbstractRaster{T,N}; kw...) where {T,N}
    _def_dim_var!(ds, A)
    attrib = _attribdict(metadata(A))
    # Set _FillValue
    eltyp = Missings.nonmissingtype(T)
    eltyp <: NCDAllowedType || throw(ArgumentError("$eltyp cannot be written to NetCDF, convert to one of $(Base.uniontypes(NCDAllowedType))"))
    if ismissing(missingval(A))
        fillval = if haskey(attrib, "_FillValue") && attrib["_FillValue"] isa eltyp
            attrib["_FillValue"]
        else
            NCD.fillvalue(eltyp)
        end
        attrib["_FillValue"] = fillval
        A = replace_missing(A, fillval)
    elseif missingval(A) isa T
        attrib["_FillValue"] = missingval(A)
    else
        missingval(A) isa Nothing || @warn "`missingval` $(missingval(A)) is not the same type as your data $T."
    end

    key = if string(name(A)) == ""
        UNNAMED_NCD_FILE_KEY
    else
        string(name(A))
    end

    dimnames = lowercase.(string.(map(name, dims(A))))
    var = NCD.defVar(ds, key, eltyp, dimnames; attrib=attrib, kw...)

    # NCDatasets needs Colon indices to write without allocations
    # TODO do this with DiskArrays broadcast ??
    var[map(_ -> Colon(), axes(A))...] = parent(read(A))

    return nothing
end

_def_dim_var!(ds::AbstractDataset, A) = map(d -> _def_dim_var!(ds, d), dims(A))
function _def_dim_var!(ds::AbstractDataset, dim::Dimension)
    dimkey = lowercase(string(name(dim)))
    haskey(ds.dim, dimkey) && return nothing
    NCD.defDim(ds, dimkey, length(dim))
    lookup(dim) isa NoLookup && return nothing

    # Shift index before conversion to Mapped
    dim = _ncdshiftlocus(dim)
    if dim isa Y || dim isa X
        dim = convertlookup(Mapped, dim)
    end
    # Attributes
    attrib = _attribdict(metadata(dim))
    _ncd_set_axis_attrib!(attrib, dim)
    # Bounds variables
    if span(dim) isa Explicit
        bounds = val(span(dim))
        boundskey = get(metadata(dim), :bounds, string(dimkey, "_bnds"))
        push!(attrib, "bounds" => boundskey)
        NCD.defVar(ds, boundskey, bounds, ("bnds", dimkey))
    end
    NCD.defVar(ds, dimkey, Vector(index(dim)), (dimkey,); attrib=attrib)
    return nothing
end

# Add axis and standard name attributes to dimension variabls
# We need to get better at guaranteeing if X/Y is actually measured in `longitude/latitude`
# CF standards requires that we specify "units" if we use these standard names
_ncd_set_axis_attrib!(atr, dim::X) = atr["axis"] = "X" # at["standard_name"] = "longitude";
_ncd_set_axis_attrib!(atr, dim::Y) = atr["axis"] = "Y" # at["standard_name"] = "latitude"; 
_ncd_set_axis_attrib!(atr, dim::Z) = (atr["axis"] = "Z"; atr["standard_name"] = "depth")
_ncd_set_axis_attrib!(atr, dim::Ti) = (atr["axis"] = "T"; atr["standard_name"] = "time")
_ncd_set_axis_attrib!(atr, dim) = nothing

_ncdshiftlocus(dim::Dimension) = _ncdshiftlocus(lookup(dim), dim)
_ncdshiftlocus(::LookupArray, dim::Dimension) = dim
function _ncdshiftlocus(lookup::AbstractSampled, dim::Dimension)
    if span(lookup) isa Regular && sampling(lookup) isa Intervals
        # We cant easily shift a DateTime value
        if eltype(dim) isa Dates.AbstractDateTime
            if !(locus(dim) isa Center)
                @warn "To save to netcdf, DateTime values should be the interval Center, rather than the $(nameof(typeof(locus(dim))))"
            end
            dim
        else
            shiftlocus(Center(), dim)
        end
    else
        dim
    end
end

_unuseddimerror(dimname) = error("Netcdf contains unused dimension $dimname")

function _ncd_eachchunk(var)
    # chunklookup, chunkvec = NCDatasets.chunking(var)
    # chunksize = chunklookup == :chunked ? Tuple(chunkvec) :
    chunksize = size(var)
    DA.GridChunks(var, chunksize)
end

function _ncd_haschunks(var)
    # chunklookup, _ = NCDatasets.chunking(var)
    # chunklookup == :chunked ? DA.Chunked() :
    DA.Unchunked()
end

# precompilation

const _NCDVar = NCDatasets.CFVariable{Union{Missing, Float32}, 3, NCDatasets.Variable{Float32, 3, NCDatasets.NCDataset}, NCDatasets.Attributes{NCDatasets.NCDataset{Nothing}}, NamedTuple{(:fillvalue, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor), Tuple{Float32, Nothing, Nothing, Nothing, Nothing, Nothing}}}

function _precompile(::Type{NCDsource})
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing

    precompile(Rasters.FileArray, (_NCDVar, String))
    precompile(layerkeys, (NCDatasets.NCDataset{Nothing},))
    precompile(dims, (_NCDVar,Symbol))
    precompile(dims, (_NCDVar,Symbol,Nothing,Nothing))
    precompile(dims, (_NCDVar,Symbol,Nothing,EPSG))
    precompile(dims, (_NCDVar,Symbol,EPSG,EPSG))
    precompile(_firstkey, (NCDatasets.NCDataset{Nothing},))
    precompile(_ncddim, (NCDatasets.NCDataset{Nothing}, Symbol, Nothing, Nothing))
    precompile(_ncddim, (NCDatasets.NCDataset{Nothing}, Symbol, Nothing, EPSG))
    precompile(_ncddim, (NCDatasets.NCDataset{Nothing}, Symbol, EPSG, EPSG))
    precompile(Raster, (NCDatasets.NCDataset{Nothing}, String, Nothing))
    precompile(Raster, (NCDatasets.NCDataset{Nothing}, String, Symbol))
    precompile(Raster, (_NCDVar, String, Symbol))

    precompile(Raster, (String,))
end

_precompile(NCDsource)

