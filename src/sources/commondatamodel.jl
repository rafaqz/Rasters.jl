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

Base.getindex(os::OpenStack{<:CDMsource}, name::Symbol) = CFDiskArray(dataset(os)[name])

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

function FileArray{source}(var::AbstractVariable, filename::AbstractString; kw...) where source<:CDMsource
    eachchunk = DA.eachchunk(var)
    haschunks = DA.haschunks(var)
    T = eltype(var)
    N = ndims(var)
    FileArray{source,T,N}(filename, size(var); eachchunk, haschunks, kw...)
end

function FileStack{source}(
    ds::AbstractDataset, filename::AbstractString;
    write::Bool=false, 
    group=nokw,
    name::NTuple{N,Symbol}, 
    vars
) where {source<:CDMsource,N}
    T = NamedTuple{name,Tuple{map(var -> Union{Missing,eltype(var)}, vars)...}}
    layersizes = map(size, vars)
    eachchunk = map(_get_eachchunk, vars)
    haschunks = map(_get_haschunks, vars)
    group = isnokw(group) ? nothing : group
    return FileStack{source,name,T}(filename, layersizes, group, eachchunk, haschunks, write)
end

function Base.open(f::Function, A::FileArray{source}; write=A.write, kw...) where source<:CDMsource
    _open(source(), filename(A); name=name(A), group=A.group, write, kw...) do var
        f(var)
    end
end

function _open(f, ::CDMsource, ds::AbstractDataset; name=nokw, group=nothing, kw...)
    g = _getgroup(ds, group)
    x = isnokw(name) ? g : CFDiskArray(g[_firstname(g, name)])
    cleanreturn(f(x))
end
_open(f, ::CDMsource, var::CFDiskArray; kw...) = cleanreturn(f(var))

# This allows arbitrary group nesting
_getgroup(ds, ::Union{Nothing,NoKW}) = ds
_getgroup(ds, group::Union{Symbol,AbstractString}) = ds.group[String(group)]
_getgroup(ds, group::Pair) = _getgroup(ds.group[String(group[1])], group[2])

function create(filename, source::CDMsource, T::Type, dims::DimTuple;
    name=nokw,
    missingval=nokw,
    metadata=nokw,
    lazy=true,
    verbose=true,
    chunks=nokw,
)
    # Create layers of zero arrays
    A = FillArrays.Zeros{T}(map(length, dims))
    rast = Raster(A, dims; name, missingval, metadata)
    write(filename, source, rast; chunks)
    return Raster(filename; metadata, source, lazy)
end

filekey(ds::AbstractDataset, name) = _firstname(ds, name)
missingval(var::AbstractDataset) = missing
missingval(var::AbstractVariable{T}) where T = missing isa T ? missing : nothing
cleanreturn(A::AbstractVariable) = Array(A)
haslayers(::CDMsource) = true
defaultcrs(::CDMsource) = EPSG(4326)
defaultmappedcrs(::CDMsource) = EPSG(4326)

function _nondimnames(ds)
    dimnames = CDM.dimnames(ds)
    toremove = if "bnds" in dimnames
        dimnames = setdiff(dimnames, ("bnds",))
        boundsnames = String[]
        for k in dimnames
            var = ds[k]
            attr = CDM.attribs(var)
            if haskey(attr, "bounds")
                push!(boundsnames, attr["bounds"])
            end
        end
        union(dimnames, boundsnames)::Vector{String}
    else
        dimnames::Vector{String}
    end
    nondim = setdiff(keys(ds), toremove)
    return nondim
end

function _layers(ds::AbstractDataset, ::NoKW=nokw, ::NoKW=nokw)
    nondim = _nondimnames(ds)
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
        names=nondim[bitinds],
        vars=vars[bitinds],
        attrs=attrs[bitinds],
    )
end
function _layers(ds::AbstractDataset, names, ::NoKW)
    vars = map(k -> ds[k], names)
    attrs = map(CDM.attribs, vars)
    (; names, vars, attrs)
end
function _layers(ds::AbstractDataset, names, group)
    _layers(ds.group[group], names, nokw)
end

function _dims(var::AbstractVariable{<:Any,N}, crs=nokw, mappedcrs=nokw) where N
    dimnames = CDM.dimnames(var)
    ntuple(Val(N)) do i
        _cdmdim(CDM.dataset(var), dimnames[i], crs, mappedcrs)
    end
end
_metadata(var::AbstractVariable; attr=CDM.attribs(var)) =
    _metadatadict(_sourcetrait(var), attr)

function _dimdict(ds::AbstractDataset, crs=nokw, mappedcrs=nokw)
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
_firstname(ds::AbstractDataset, name) = Symbol(name)
function _firstname(ds::AbstractDataset, name::NoKW=nokw)
    names = _nondimnames(ds)
    if length(names) > 0
        Symbol(first(names))
    else
        throw(ArgumentError("No non-dimension layers found in dataset with keys: $(keys(ds))"))
    end
end

function _cdmdim(ds, dimname::Key, crs=nokw, mappedcrs=nokw)
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
    for name in keys(ds)
        var = ds[name]
        dimnames = CDM.dimnames(var)
        if dimname in dimnames
            return size(var)[findfirst(==(dimname), dimnames)]
        end
    end
    return nothsng
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
    return DD.basetypeof(DD.name2dim(Symbol(dimname)))
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
    # Detect lat/lon
    span, sampling = if eltype(index) <: Union{Missing,Dates.AbstractTime}
        _cdmperiod(index, metadata)
    else
        _cdmspan(index, order)
    end
    # We only use Explicit if the span is not Regular
    # This is important for things like rasterizatin and conversion 
    # to gdal to be easy, and selectors are faster.
    # TODO are there any possible floating point errors from this?
    if haskey(CDM.attribs(var), "bounds")
        span, sampling = if isregular(span)
            span, Intervals(Center())
        else
            boundskey = var.attrib["bounds"]
            boundsmatrix = Array(ds[boundskey])
            Explicit(boundsmatrix), Intervals(Center())
        end
    end

    # We cant yet check CF standards crs, but we can at least check for units in lat/lon 
    if isnokw(mappedcrs) && get(metadata, "units", "") in ("degrees_north", "degrees_east")
        mappedcrs = EPSG(4326)
    end
    # Additionally, crs and mappedcrs should be identical for Regular lookups
    if isnokw(crs) && span isa Regular
        crs = mappedcrs
    end
    return _cdmlookup(D, index, order, span, sampling, metadata, crs, mappedcrs)
end
# For X and Y use a Mapped <: AbstractSampled lookup
function _cdmlookup(
    D::Type{<:Union{<:XDim,<:YDim}}, index, order::Order, span, sampling, metadata, crs, mappedcrs
)
    # If the index is regularly spaced and there is no crs
    # then there is probably just one crs - the mappedcrs
    crs = if isnokw(crs) && span isa Regular
        mappedcrs
    else
        crs
    end
    dim = DD.basetypeof(D)()
    return Mapped(index; order, span, sampling, metadata, crs, mappedcrs, dim)
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
    length(index) == 1 && return Regular(zero(eltype(index))), Points()

    step = if eltype(index) <: AbstractFloat
        # Calculate step, avoiding as many floating point errors as possible
        st = Base.step(Base.range(Float64(first(index)), Float64(last(index)); length = length(index)))
        st_rd = round(st, digits = Base.floor(Int,-log10(eps(eltype(index))))) # round to nearest digit within machine epsilon
        st_rd ≈ st ? st_rd : st # keep the rounded number if it is very close to the original
    else
        index[2] - index[1]
    end
    
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
            return Irregular(bounds), Points()
        end
    end
    # Otherwise regular
    return Regular(step), Points()
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
            mo = Month(vals[2])
            d = Day(vals[3])
            h = Hour(vals[4])
            mi = Minute(vals[5])
            s = Second(vals[6])
            compound = sum((y, mo, d, h, mi, s))
            if length(compound.periods) == 0
                return nothing
            elseif length(compound.periods) == 1
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
