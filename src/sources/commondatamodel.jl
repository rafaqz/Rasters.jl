const CDM = CommonDataModel

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

struct CFVariable{T,N,TV,TA,TSA} <: CDM.AbstractVariable{T,N}
    var::CDM.CFVariable{T,N,TV,TA,TSA}
end
@implement_diskarray CFVariable

Base.parent(A::CFVariable) = A.var

function DiskArrays.readblock!(A::CFVariable, aout, i::AbstractUnitRange...)
    aout[i...] = getindex(parent(A), i...)
end

function DiskArrays.writeblock!(A::CFVariable, data, i::AbstractUnitRange...)
    setindex!(parent(A), data, i...)
    return data
end
_open(f, ::Type{<:CDMsource}, var::CFVariable; kw...) = cleanreturn(f(var))

# Methods from CommonDataModel
for method in (:size, :name, :dimnames, :dataset, :attribnames)
    @eval begin
        CDM.$(method)(var::CFVariable) = CDM.$(method)(parent(var))
    end
end

for method in (:attrib, :dim)
    @eval begin
        CDM.$(method)(var::CFVariable, name::CDM.SymbolOrString) = CDM.$(method)(parent(var), name)
    end
end

Base.getindex(os::OpenStack{<:CDMsource}, key::Symbol) = CFVariable(dataset(os)[key])

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
                return Raster(CFVariable(ds[key]), filename, key; source=CDMsource, kw...)
            end
        end
        throw(ArgumentError("dataset at $filename has no array variables"))
    else
       return Raster(CFVariable(ds[key]), filename, key; kw...)
    end
end

_firstkey(ds::AbstractDataset, key::Nothing=nothing) = Symbol(first(layerkeys(ds)))
_firstkey(ds::AbstractDataset, key) = Symbol(key)

function FileArray(var::AbstractVariable, filename::AbstractString; kw...)
    source = _sourcetype(filename)
    da = RasterDiskArray{source}(var)
    size_ = size(da)
    eachchunk = DA.eachchunk(da)
    haschunks = DA.haschunks(da)
    T = eltype(var)
    N = length(size_)
    FileArray{source,T,N}(filename, size_; eachchunk, haschunks, kw...)
end

function Base.open(f::Function, A::FileArray{source}; write=A.write, kw...) where source <: CDMsource
    _open(source, filename(A); key=key(A), write, kw...) do var
        f(var)
    end
end

function create(filename, source::Type{<:CDMsource}, T::Union{Type,Tuple}, dims::DimTuple;
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
    write(filename, source, Raster(first(layers)))
    return Raster(filename; source=source, lazy)
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
DD.metadata(var::AbstractVariable) = _metadatadict(CDMsource, _attrib(var))

function DD.layerdims(ds::AbstractDataset)
    keys = Tuple(layerkeys(ds))
    dimtypes = map(keys) do key
        DD.layerdims(ds[string(key)])
    end
    NamedTuple{map(Symbol, keys)}(dimtypes)
end
function DD.layerdims(var::AbstractVariable)
    map(CDM.dimnames(var)) do dimname
        _ncddim(_dataset(var), dimname)
    end |> Tuple    
end

function DD.layermetadata(ds::AbstractDataset)
    keys = Tuple(layerkeys(ds))
    dimtypes = map(keys) do k
        var = ds[k]
        md = DD.metadata(var)
        if haskey(_attrib(var), "grid_mapping")
            md["grid_mapping"] = Dict(_attrib(ds[_attrib(var)["grid_mapping"]]))
        end
        md
    end
    NamedTuple{map(Symbol, keys)}(dimtypes)
end

missingval(var::AbstractVariable{T}) where T = missing isa T ? missing : nothing

function layerkeys(ds::AbstractDataset)
    dimkeys = _dimkeys(ds)
    toremove = if "bnds" in dimkeys
        dimkeys = setdiff(dimkeys, ("bnds",))
        boundskeys = String[]
        for k in dimkeys
            var = ds[k]
            if haskey(_attrib(var), "bounds")
                push!(boundskeys, _attrib(var)["bounds"])
            end
        end
        union(dimkeys, boundskeys)::Vector{String}
    else
        dimkeys::Vector{String}
    end
    nondim = setdiff(keys(ds), toremove)
    grid_mapping = String[]
    for k in nondim
        var = ds[k]
        if haskey(_attrib(var), "grid_mapping")
            push!(grid_mapping, _attrib(var)["grid_mapping"])
        end
    end
    nondim = setdiff(nondim, grid_mapping)
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

function _open(f, ::Type{<:CDMsource}, ds::AbstractDataset; key=nothing, kw...)
    x = key isa Nothing ? ds : CFVariable(ds[_firstkey(ds, key)])
    cleanreturn(f(x))
end
# _open(f, ::Type{<:CDMsource}, var::CDM.CFVariable; kw...) = cleanreturn(f(CFVariable(var)))

# Utils ########################################################################

cleanreturn(A::CFVariable) = Array(A)

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


