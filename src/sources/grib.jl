const GRIB_DIM_MAP = Dict(
    "latitude" => Y,
    "longitude" => X,
    "valid_time" => Ti,
    "level" => Z,
    "vertical" => Z,
    "x" => X,
    "y" => Y,
    "z" => Z,
)

haslayers(::Type{GRIBfile}) = true
defaultcrs(::Type{GRIBfile}) = EPSG(4326)
defaultmappedcrs(::Type{GRIBfile}) = EPSG(4326)

# Implement some methods on Dataset and Variable that are needed for conversion to Raster
Base.Array(var::GDS.Variable) = Array(var.values)
DA.readblock!(var::GDS.Variable, args...) = DA.readblock!(var.values, args...)

# Raster ########################################################################

function Raster(ds::GDS.Dataset, filename::AbstractString, key=nothing; kw...)
    key = _firstkey(ds, key)
    Raster(ds[key], filename, key; kw...)
end

_firstkey(ds::GDS.Dataset, key::Nothing=nothing) = Symbol(first(layerkeys(ds)))
_firstkey(ds::GDS.Dataset, key) = Symbol(key)

_dimkeys(ds::GDS.Dataset) = collect(keys(ds.dims))
_dimnames(var::GDS.Variable) = collect(keys(var.dims))


function FileArray(var::GDS.Variable, filename::AbstractString; kw...)
    da = RasterDiskArray{GRIBfile}(var.values)
    size_ = size(da)
    eachchunk = DA.eachchunk(da)
    haschunks = DA.haschunks(da)
    T = eltype(var)
    N = length(size_)
    FileArray{GRIBfile,T,N}(filename, size_; eachchunk, haschunks, kw...)
end


# DimensionalData methods for CfGRIB.Dataset types ###############################

function DD.dims(ds::GDS.Dataset, crs=nothing, mappedcrs=nothing)
    map(_dimkeys(ds)) do key
        _dsdim(ds, key, crs, mappedcrs)
    end |> Tuple
end
function DD.dims(var::GDS.Variable, crs=nothing, mappedcrs=nothing)
    names = _dimnames(var)
    ds = var.ds
    map(names) do name
        _dsdim(ds, name, crs, mappedcrs)
    end |> Tuple
end

DD.metadata(ds::GDS.Dataset) = _metadatadict(GRIBfile, ds.attrib)
DD.metadata(var::GDS.Variable) = _metadatadict(GRIBfile, var.attrib)

function DD.layerdims(ds::GDS.Dataset)
    keys = Tuple(layerkeys(ds))
    dimtypes = map(keys) do key
        DD.layerdims(ds[string(key)])
    end
    NamedTuple{map(Symbol, keys)}(dimtypes)
end
function DD.layerdims(var::GDS.Variable)
    map(_dimnames(var)) do dimname
        _gdsdimtype(dimname)()
    end |> Tuple
end

DD.layermetadata(ds::GDS.Dataset) = _layermetadata(ds, Tuple(layerkeys(ds)))
function _layermetadata(ds::GDS.Dataset, keys)
    dimtypes = map(k -> DD.metadata(ds[string(k)]), keys)
    NamedTuple{map(Symbol, keys)}(dimtypes)
end

missingval(var::GDS.Variable) = missing
missingval(ds::GDS.Dataset) = missing

function layerkeys(ds::GDS.Dataset)
    # dimension_keys = [k for (k, v) in ds.variables if length(size(v.data)) <= 1]
    GDS.getlayersname(ds)
end

function FileStack{GRIBfile}(ds::GDS.Dataset, filename::AbstractString; write=false, keys)
    keys = map(Symbol, keys isa Nothing ? collect(layerkeys(ds)) : keys) |> Tuple
    type_size_ec_hc = map(keys) do key
        var = ds[string(key)]
        Union{Missing,eltype(var)}, size(var), _gds_eachchunk(var), _gds_haschunks(var)
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
    # CFG.Dataset(filename) do ds
    #     _open(f, GRIBfile, ds; kw...)
    # end
    ds = GDS.Dataset(filename)
    _open(f, GRIBfile, ds; kw...)
end
function _open(f, ::Type{GRIBfile}, ds::GDS.Dataset; key=nothing, kw...)
    cleanreturn(f(key isa Nothing ? ds : ds[_firstkey(ds, key)]))
end
_open(f, ::Type{GRIBfile}, var::GDS.Variable; kw...) = cleanreturn(f(var))

cleanreturn(A::GDS.Variable) = Array(A.values)

function _dsdim(ds::GDS.Dataset, dimname::Key, crs=nothing, mappedcrs=nothing)
    D = _gdsdimtype(dimname)
    lookup = _dslookup(ds, dimname, D, crs, mappedcrs)
    return D(lookup)
end

function _dslookup(ds::GDS.Dataset, dimname, D::Type, crs, mappedcrs)
    dvar = ds[dimname]
    index = dvar[:]
    metadata = _metadatadict(GRIBfile, dvar.attrib)
    return _dslookup(ds, dimname, D, index, metadata, crs, mappedcrs)
end

function _dslookup(
    ds::GDS.Dataset, dimname, D, index::AbstractArray{<:Union{Number,Dates.AbstractTime}},
    metadata, crs, mappedcrs
)
    # Assume the locus is at the center of the cell if boundaries aren't provided.
    # http://cfconventions.org/cf-conventions/cf-conventions.html#cell-boundaries
    order = LA.orderof(index)
    # var = NCD.variable(ds, dimname)
    if dimname in ["time", "valid_time"]
        # We consider the epoch 1970-01-01T00:00:00, as it appears to be in gribs files
        # dates = Dates.unix2datetime.(index)
        dates = Second.(index) .+ GDS.DEFAULT_EPOCH
        steps = unique(dates[2:end] .- dates[1:end-1])
        if length(steps) == 1
            span, sampling = Regular(steps[1]), Points()
        else
            span, sampling = Irregular((extrema(dates)...)), Points()
        end
        index = dates
        return _dslookup(D, index, order, span, sampling, metadata, crs, mappedcrs)
    else
        span, sampling = _dsspan(index, order), Points()
        return _dslookup(D, index, order, span, sampling, metadata, crs, mappedcrs)
    end
end

# Find the matching dimension constructor. If its an unknown name
# use the generic Dim with the dim name as type parameter
_gdsdimtype(dimname) = haskey(GRIB_DIM_MAP, dimname) ? GRIB_DIM_MAP[dimname] : DD.basetypeof(DD.key2dim(Symbol(dimname)))

_gds_eachchunk(var) = DA.eachchunk(var.values)
_gds_haschunks(var) = DA.haschunks(var.values)