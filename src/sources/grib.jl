const CFG = CfGRIB

const GRIB_DIMMAP = Dict(
    "lat" => Y,
    "latitude" => Y,
    "lon" => X,
    "long" => X,
    "longitude" => X,
    "valid_time" => Ti,
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

haslayers(::Type{GRIBfile}) = true
defaultcrs(::Type{GRIBfile}) = EPSG(4326)
defaultmappedcrs(::Type{GRIBfile}) = EPSG(4326)

# Implement some methods on DataSet and Variable that are needed for conversion to Raster

Base.getindex(ds::CFG.DataSet, index) = getindex(ds.variables, String(index))
Base.getindex(var::CFG.Variable, index...) = getindex(var.data, index...)
# Base.getindex(var::CFG.Variable, index) = getindex(var, String(index))
Base.haskey(ds::CFG.DataSet, key) = haskey(ds.variables, String(key))
Base.Array(var::CFG.Variable) = Array(var.data)
Base.size(var::CFG.Variable) = size(var.data)
DA.readblock!(var::CFG.Variable, args...) = DA.readblock!(var.data, args...)

# Raster ########################################################################

function Raster(ds::CFG.DataSet, filename::AbstractString, key=nothing; kw...)
    key = _firstkey(ds, key)
    Raster(ds[key], filename, key; kw...)
end

_firstkey(ds::CFG.DataSet, key::Nothing=nothing) = Symbol(first(layerkeys(ds)))
_firstkey(ds::CFG.DataSet, key) = Symbol(key)

_dimkeys(ds::CFG.DataSet) = collect(keys(ds.dimensions))
_dimnames(var::CFG.Variable) = var.dimensions


function FileArray(var::CFG.Variable, filename::AbstractString; kw...)
    da = RasterDiskArray{GRIBfile}(var.data)
    size_ = size(da)
    eachchunk = DA.eachchunk(da)
    haschunks = DA.haschunks(da)
    T = eltype(var)
    N = length(size_)
    FileArray{GRIBfile,T,N}(filename, size_; eachchunk, haschunks, kw...)
end


# DimensionalData methods for CfGRIB.DataSet types ###############################

function DD.dims(ds::CFG.DataSet, crs=nothing, mappedcrs=nothing)
    map(_dimkeys(ds)) do key
        _cfgdim(ds, key, crs, mappedcrs)
    end |> Tuple
end
function DD.dims(var::CFG.Variable, crs=nothing, mappedcrs=nothing)
    names = _dimnames(var)
    map(names) do name
        # Not good since we have to rebuild the DataSet each time! Should modify CfGRIB to keep a reference to the DataSet in the Variable
        ds = CFG.DataSet(var.data.grib_path)
        _cfgdim(ds, name, crs, mappedcrs)
    end |> Tuple
end

DD.metadata(ds::CFG.DataSet) = Metadata{GRIBfile}(LA.metadatadict(ds.attributes))
DD.metadata(var::CFG.Variable) = Metadata{GRIBfile}(LA.metadatadict(var.attributes))

function DD.layerdims(ds::CFG.DataSet)
    keys = Tuple(layerkeys(ds))
    dimtypes = map(keys) do key
        DD.layerdims(ds[string(key)])
    end
    NamedTuple{map(Symbol, keys)}(dimtypes)
end
function DD.layerdims(var::CFG.Variable)
    map(_dimnames(var)) do dimname
        _cfgdimtype(dimname)()
    end
end

DD.layermetadata(ds::CFG.DataSet) = _layermetadata(ds, Tuple(layerkeys(ds)))
function _layermetadata(ds::CFG.DataSet, keys)
    dimtypes = map(k -> DD.metadata(ds[string(k)]), keys)
    NamedTuple{map(Symbol, keys)}(dimtypes)
end

missingval(var::CFG.Variable) = missing
missingval(ds::CFG.DataSet) = missing

function layerkeys(ds::CFG.DataSet)
    dimension_keys = [k for (k, v) in ds.variables if length(size(v.data)) <= 1]
    return setdiff(keys(ds.variables), dimension_keys)
end

function FileStack{GRIBfile}(ds::CFG.DataSet, filename::AbstractString; write=false, keys)
    keys = map(Symbol, keys isa Nothing ? collect(layerkeys(ds)) : keys) |> Tuple
    type_size_ec_hc = map(keys) do key
        var = ds[string(key)]
        Union{Missing,eltype(var)}, size(var), _cfg_eachchunk(var), _cfg_haschunks(var)
    end
    layertypes = map(x->x[1], type_size_ec_hc)
    layersizes = map(x->x[2], type_size_ec_hc)
    eachchunk = map(x->x[3], type_size_ec_hc)
    haschunks = map(x->x[4], type_size_ec_hc)
    return FileStack{GRIBfile,keys}(filename, layertypes, layersizes, eachchunk, haschunks, write)
end

# Utils ########################################################################

function _open(f, ::Type{GRIBfile}, filename::AbstractString; write=false, kw...)
    isfile(filename) || _filenotfound_error(filename)
    mode = write ? "a" : "r"
    # CFG.DataSet(filename) do ds
    #     _open(f, GRIBfile, ds; kw...)
    # end
    ds = CFG.DataSet(filename)
    _open(f, GRIBfile, ds; kw...)
end
function _open(f, ::Type{GRIBfile}, ds::CFG.DataSet; key=nothing, kw...)
    cleanreturn(f(key isa Nothing ? ds : ds[_firstkey(ds, key)]))
end
_open(f, ::Type{GRIBfile}, var::CFG.Variable; kw...) = cleanreturn(f(var))

cleanreturn(A::CFG.Variable) = Array(A.data)

function _cfgdim(ds, dimname::Key, crs=nothing, mappedcrs=nothing)
    if haskey(ds, dimname)
        D = _cfgdimtype(dimname)
        lookup = _cfglookup(ds, dimname, D, crs, mappedcrs)
        return D(lookup)
    else
        # The var doesn't exist. Maybe its `complex` or some other marker,
        # so make it a custom `Dim` with `NoLookup`
        # len = _cfgfinddimlen(ds, dimname)
        # len === nothing && _unuseddimerror(dimname)
        # return Dim{Symbol(dimname)}(NoLookup(Base.OneTo(len)))
        error("GRIB contains unused dim")
    end
end

# Generate a `LookupArray` from a netcdf dim.
function _cfglookup(ds::CFG.DataSet, dimname, D, crs, mappedcrs)
    dvar = ds[dimname]
    index = dvar[:]
    metadata = Metadata{GRIBfile}(LA.metadatadict(dvar.attributes))
    return _cfglookup(ds, dimname, D, index, metadata, crs, mappedcrs)
end
# For unknown types we just make a Categorical lookup
function _cfglookup(ds::CFG.DataSet, dimname, D, index::AbstractArray, metadata, crs, mappedcrs)
    Categorical(index; metadata=metadata)
end
# For Number and AbstractTime we generate order/span/sampling
function _cfglookup(
    ds::CFG.DataSet, dimname, D, index::AbstractArray{<:Union{Number,Dates.AbstractTime}},
    metadata, crs, mappedcrs
)
    # Assume the locus is at the center of the cell if boundaries aren't provided.
    # http://cfconventions.org/cf-conventions/cf-conventions.html#cell-boundaries
    order = _cfgorder(index)
    # var = NCD.variable(ds, dimname)
    if dimname in ["time", "valid_time"]
        # We consider the epoch 1970-01-01T00:00:00, as it appears to be in gribs files
        dates = Dates.unix2datetime.(index)
        steps = unique(dates[2:end] .- dates[1:end-1])
        if length(steps) == 1
            span, sampling = Regular(steps[1]), Points()
        else
            span, sampling = Irregular((extrema(dates)...)), Points()
        end
        return _cfglookup(D, index, order, span, sampling, metadata, crs, mappedcrs)
    # elseif eltype(index) <: Dates.AbstractTime
    #     span, sampling = _cfgperiod(index, metadata)
    #     return _cfglookup(D, index, order, span, sampling, metadata, crs, mappedcrs)
    else
        # Use the same function as for NetCDF files
        span, sampling = _ncdspan(index, order), Points()
        return _cfglookup(D, index, order, span, sampling, metadata, crs, mappedcrs)
    end
end
# For X and Y use a Mapped <: AbstractSampled lookup
function _cfglookup(
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
function _cfglookup(D::Type{<:Band}, index, order::Order, span, sampling, metadata, crs, mappedcrs)
    Categorical(index, order, metadata)
end
# Otherwise use a regular Sampled lookup
function _cfglookup(D::Type, index, order::Order, span, sampling, metadata, crs, mappedcrs)
    Sampled(index, order, span, sampling, metadata)
end

function _cfgorder(index)
    index[end] > index[1] ? ForwardOrdered() : ReverseOrdered()
end

# Find the matching dimension constructor. If its an unknown name
# use the generic Dim with the dim name as type parameter
_cfgdimtype(dimname) = haskey(GRIB_DIMMAP, dimname) ? GRIB_DIMMAP[dimname] : DD.basetypeof(DD.key2dim(Symbol(dimname)))

_cfg_eachchunk(var) = DA.eachchunk(var.data)
_cfg_haschunks(var) = DA.haschunks(var.data)