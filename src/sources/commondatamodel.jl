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


# CFDiskArray ########################################################################

struct CFDiskArray{T,N,TV,TA,TSA} <: DiskArrays.AbstractDiskArray{T,N}
    var::CDM.CFVariable{T,N,TV,TA,TSA}
end

# Rasters methods
FileArray{source}(var::CFDiskArray, filename::AbstractString; kw...) where source =
    FileArray{source}(parent(var), filename; kw...)

cleanreturn(A::CFDiskArray) = Array(A)
missingval(A::CFDiskArray) = missingval(parent(A))

# DimensionalData methods
_dims(var::CFDiskArray, args...) = _dims(parent(var), args...)
_metadata(var::CFDiskArray, args...) = _metadata(parent(var), args...)

# Base methods
Base.parent(A::CFDiskArray) = A.var

Base.getindex(os::OpenStack{<:CDMsource}, key::Symbol) = CFDiskArray(dataset(os)[key])

# DiskArrays.jl methods
function DiskArrays.readblock!(A::CFDiskArray, aout, i::AbstractUnitRange...)
    aout .= getindex(parent(A), i...)
end
function DiskArrays.writeblock!(A::CFDiskArray, data, i::AbstractUnitRange...)
    setindex!(parent(A), data, i...)
    return data
end

# We have to dig down to find the chunks as they are not immplemented
# in the CDM, but they are in their internal objects.
DiskArrays.eachchunk(var::CFDiskArray) = _get_eachchunk(var)
DiskArrays.haschunks(var::CFDiskArray) = _get_haschunks(var)

_get_eachchunk(var::CFDiskArray) = _get_eachchunk(parent(var))
_get_eachchunk(var::CDM.CFVariable) = _get_eachchunk(var.var)
_get_haschunks(var::CFDiskArray) = _get_haschunks(parent(var))
_get_haschunks(var::CDM.CFVariable) = _get_haschunks(var.var)

_sourcetrait(var::CFDiskArray) = _sourcetrait(parent(var))
_sourcetrait(var::CDM.CFVariable) = _sourcetrait(var.var)

# CommonDataModel.jl methods
for method in (:size, :name, :dimnames, :dataset, :attribnames)
    @eval begin
        CDM.$(method)(var::CFDiskArray) = CDM.$(method)(parent(var))
    end
end

for method in (:attrib, :dim)
    @eval begin
        CDM.$(method)(var::CFDiskArray, name::CDM.SymbolOrString) = CDM.$(method)(parent(var), name)
    end
end

# Rasters methods for CDM types ###############################

# This is usually called inside a closure and cleaned up in `cleanreturn`
function Raster(ds::AbstractDataset, filename::AbstractString, key::Nothing=nothing;
    source=nothing, kw...
)
    source = isnothing(source) ? _sourcetrait(filename) : _sourcetrait(source)
    # Find the first valid variable
    layers = _layers(ds)
    for (key, var) in zip(layers.keys, layers.vars)
        if ndims(var) > 0
            @info "No `name` or `key` keyword provided, using first valid layer with name `:$key`"
            return Raster(CFDiskArray(var), filename, key; source, kw...)
        end
    end
    throw(ArgumentError("dataset at $filename has no array variables"))
end
function Raster(ds::AbstractDataset, filename::AbstractString, key::Union{AbstractString,Symbol};
    source=nothing, kw...
)
    return Raster(CFDiskArray(ds[key]), filename, key; source)
end

function FileArray{source}(var::AbstractVariable, filename::AbstractString; kw...) where source<:CDMsource
    eachchunk = DA.eachchunk(var)
    haschunks = DA.haschunks(var)
    T = eltype(var)
    N = ndims(var)
    FileArray{source,T,N}(filename, size(var); eachchunk, haschunks, kw...)
end

function FileStack{source}(
    ds::AbstractDataset, filename::AbstractString;
    write::Bool=false, keys::NTuple{N,Symbol}, vars
) where {source<:CDMsource,N}
    layertypes = map(var -> Union{Missing,eltype(var)}, vars)
    layersizes = map(size, vars)
    eachchunk = map(_get_eachchunk, vars)
    haschunks = map(_get_haschunks, vars)
    return FileStack{source,keys}(filename, layertypes, layersizes, eachchunk, haschunks, write)
end

function Base.open(f::Function, A::FileArray{source}; write=A.write, kw...) where source<:CDMsource
    _open(source(), filename(A); key=key(A), write, kw...) do var
        f(var)
    end
end

function _open(f, ::CDMsource, ds::AbstractDataset; key=nothing, kw...)
    x = key isa Nothing ? ds : CFDiskArray(ds[_firstkey(ds, key)])
    cleanreturn(f(x))
end
_open(f, ::CDMsource, var::CFDiskArray; kw...) = cleanreturn(f(var))
# _open(f, ::CDMsource, var::CDM.CFVariable; kw...) = cleanreturn(f(CFDiskArray(var)))

# TODO fix/test this for RasterStack
function create(filename, source::CDMsource, T::Union{Type,Tuple}, dims::DimTuple;
    name=:layer1,
    keys=(name,),
    layerdims=map(_ -> dims, keys),
    missingval=nothing,
    metadata=NoMetadata(),
    lazy=true,
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

missingval(var::AbstractDataset) = missing
missingval(var::AbstractVariable{T}) where T = missing isa T ? missing : nothing
cleanreturn(A::AbstractVariable) = Array(A)
haslayers(::CDMsource) = true
defaultcrs(::CDMsource) = EPSG(4326)
defaultmappedcrs(::CDMsource) = EPSG(4326)

function _layers(ds::AbstractDataset, ::Nothing=nothing)
    dimkeys = CDM.dimnames(ds)
    toremove = if "bnds" in dimkeys
        dimkeys = setdiff(dimkeys, ("bnds",))
        boundskeys = String[]
        for k in dimkeys
            var = ds[k]
            attr = CDM.attribs(var)
            if haskey(attr, "bounds")
                push!(boundskeys, attr["bounds"])
            end
        end
        union(dimkeys, boundskeys)::Vector{String}
    else
        dimkeys::Vector{String}
    end
    nondim = setdiff(keys(ds), toremove)
    grid_mapping = String[]
    vars = map(k -> ds[k], nondim)
    attrs = map(CDM.attribs, vars)
    for attr in attrs
        if haskey(attr, "grid_mapping")
            push!(grid_mapping, attr["grid_mapping"])
        end
    end
    bitinds = map(!in(grid_mapping), nondim)
    (;
        keys=nondim[bitinds],
        vars=vars[bitinds],
        attrs=attrs[bitinds],
    )
end
function _layers(ds::AbstractDataset, keys)
    vars = map(k -> ds[k], keys)
    attrs = map(CDM.attribs, vars)
    (; keys, vars, attrs)
end

function _dims(var::AbstractVariable{<:Any,N}, crs=nothing, mappedcrs=nothing) where N
    dimnames = CDM.dimnames(var)
    ntuple(Val(N)) do i
        _cdmdim(CDM.dataset(var), dimnames[i], crs, mappedcrs)
    end
end
_metadata(var::AbstractVariable; attr=CDM.attribs(var)) =
    _metadatadict(_sourcetrait(var), attr)

function _dimdict(ds::AbstractDataset, crs=nothing, mappedcrs=nothing)
    dimdict = Dict{String,Dimension}()
    for dimname in CDM.dimnames(ds)
        dimdict[dimname] = _cdmdim(ds, dimname, crs, mappedcrs)
    end
    return dimdict
end
function _dims(ds::AbstractDataset, dimdict::Dict)
    map(CDM.dimnames(ds)) do dimname
        dimdict[dimname]
    end |> Tuple
end
_metadata(ds::AbstractDataset; attr=CDM.attribs(ds)) =
    _metadatadict(_sourcetrait(ds), attr)
function _layerdims(ds::AbstractDataset; layers, dimdict)
    map(layers.vars) do var
        map(CDM.dimnames(var)) do dimname
            basedims(dimdict[dimname])
        end |> Tuple
    end
end
function _layermetadata(ds::AbstractDataset; layers)
    map(layers.attrs) do attr
        md = _metadatadict(_sourcetrait(ds), attr)
        if haskey(attr, "grid_mapping")
            md["grid_mapping"] = Dict(CDM.attribs(ds[attr["grid_mapping"]]))
        end
        md
    end
end


# Utils ########################################################################

# TODO dont load all keys here with _layers
_firstkey(ds::AbstractDataset, key::Nothing=nothing) = Symbol(first(_layers(ds).keys))
_firstkey(ds::AbstractDataset, key) = Symbol(key)

function _cdmdim(ds, dimname::Key, crs=nothing, mappedcrs=nothing)
    if haskey(ds, dimname)
        var = ds[dimname]
        D = _cdmdimtype(CDM.attribs(var), dimname)
        lookup = _cdmlookup(ds, dimname, D, crs, mappedcrs)
        return D(lookup)
    else
        # The var doesn't exist. Maybe its `complex` or some other marker,
        # so make it a custom `Dim` with `NoLookup`
        len = _cdmfinddimlen(ds, dimname)
        len === nothing && _unuseddimerror()
        lookup = NoLookup(Base.OneTo(len))
        D = _cdmdimtype(NoMetadata(), dimname)
        return D(lookup)
    end
end

function _cdmfinddimlen(ds, dimname)
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
function _cdmdimtype(attrib, dimname)
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

# _cdmlookup
# Generate a `Lookup` from a nCDM dim.
function _cdmlookup(ds::AbstractDataset, dimname, D::Type, crs, mappedcrs)
    var = ds[dimname]
    index = var[:]
    attr = CDM.attribs(var)
    metadata = _metadatadict(_sourcetrait(ds), attr)
    return _cdmlookup(ds, var, attr, dimname, D, index, metadata, crs, mappedcrs)
end
# For unknown types we just make a Categorical lookup
function _cdmlookup(ds::AbstractDataset, var, attr, dimname, D::Type, index::AbstractArray, metadata, crs, mappedcrs)
    Categorical(index; order=Unordered(), metadata=metadata)
end
# For Number and AbstractTime we generate order/span/sampling
# We need to include `Missing` in unions in case `_FillValue` is used
# on coordinate variables in a file and propagates here.
function _cdmlookup(
    ds::AbstractDataset, var, attr, dimname,
    D::Type, index::AbstractArray{<:Union{Missing,Number,Dates.AbstractTime}},
    metadata, crs, mappedcrs
)
    # Assume the locus is at the center of the cell if boundaries aren't provided.
    # http://cfconventions.org/cf-conventions/cf-conventions.html#cell-boundaries
    order = LA.orderof(index)
    var = ds[dimname]
    if haskey(CDM.attribs(var), "bounds")
        boundskey = var.attrib["bounds"]
        boundsmatrix = Array(ds[boundskey])
        span, sampling = Explicit(boundsmatrix), Intervals(Center())
        return _cdmlookup(D, index, order, span, sampling, metadata, crs, mappedcrs)
    elseif eltype(index) <: Union{Missing,Dates.AbstractTime}
        span, sampling = _cdmperiod(index, metadata)
        return _cdmlookup(D, index, order, span, sampling, metadata, crs, mappedcrs)
    else
        span, sampling = _cdmspan(index, order), Points()
        return _cdmlookup(D, index, order, span, sampling, metadata, crs, mappedcrs)
    end
end
# For X and Y use a Mapped <: AbstractSampled lookup
function _cdmlookup(
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
function _cdmlookup(D::Type{<:Band}, index, order::Order, span, sampling, metadata, crs, mappedcrs)
    Categorical(index, order, metadata)
end
# Otherwise use a regular Sampled lookup
function _cdmlookup(D::Type, index, order::Order, span, sampling, metadata, crs, mappedcrs)
    Sampled(index, order, span, sampling, metadata)
end

function _cdmspan(index, order)
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
function _cdmperiod(index, metadata::Metadata{<:CDMsource})
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
        vals = map(x -> parse(Int, x), mtch.captures)
        if length(vals) == 6
            y = Year(vals[1])
            m = Month(vals[2])
            d = Day(vals[3])
            h = Hour(vals[4])
            m = Minute(vals[5])
            s = Second(vals[6])
            compound = sum(y, m, d, h, m, s)
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

# Add axis and standard name attributes to dimension variabls
# We need to get better at guaranteeing if X/Y is actually measured in `longitude/latitude`
# CF standards requires that we specify "units" if we use these standard names
_cdm_set_axis_attrib!(atr, dim::X) = atr["axis"] = "X" # at["standard_name"] = "longitude";
_cdm_set_axis_attrib!(atr, dim::Y) = atr["axis"] = "Y" # at["standard_name"] = "latitude";
_cdm_set_axis_attrib!(atr, dim::Z) = (atr["axis"] = "Z"; atr["standard_name"] = "depth")
_cdm_set_axis_attrib!(atr, dim::Ti) = (atr["axis"] = "T"; atr["standard_name"] = "time")
_cdm_set_axis_attrib!(atr, dim) = nothing

_cdmshiftlocus(dim::Dimension) = _cdmshiftlocus(lookup(dim), dim)
_cdmshiftlocus(::Lookup, dim::Dimension) = dim
function _cdmshiftlocus(lookup::AbstractSampled, dim::Dimension)
    # TODO improve this
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

_unuseddimerror(dimname) = error("Dataset contains unused dimension $dimname")
