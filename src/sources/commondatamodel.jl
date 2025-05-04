const CDM = CommonDataModel

const CDMallowedType = Union{Int8,UInt8,Int16,UInt16,Int32,UInt32,Int64,UInt64,Float32,Float64,Char,String}

const CF_DEGREES_NORTH = ("degrees_north", "degree_north", "degree_N", "degrees_N", "degreeN", "degreesN")
const CF_DEGREES_EAST = ("degrees_east", "degree_east", "degree_E", "degrees_E", "degreeE", "degreesE")

const UNNAMED_FILE_KEY = "unnamed"

const CDM_AXIS_MAP = Dict(
    "X" => X,
    "Y" => Y,
    "Z" => Z,
    "T" => Ti,
)

const CDM_STANDARD_NAME_MAP = Dict(
    "longitude" => X,
    "latitude" => Y,
    "height" => Z,
    "depth" => Z,
    "time" => Ti,
)

# These are mostly a hack for when
# standard_name or axis are missing
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

# `Source`` from variables and datasets
sourcetrait(var::CDM.CFVariable) = sourcetrait(var.var)
# Dataset constructor from `Source`
sourceconstructor(source::Source) = sourceconstructor(typeof(source))
# Function to check filename
checkfilename(s::CDMsource, filename) = throw(BackendException(s))
# CDM datasets are always multilayer
check_multilayer_dataset(ds::CDM.AbstractDataset) = true

# Find and check write modes
function checkwritemode(::CDMsource, filename, append::Bool, force::Bool)
    if append
        isfile(filename) ? "a" : "c"
    else
        check_can_write(filename, force)
        "c"
    end
end
# Mode to open file in - read or append
openmode(write::Bool) = write ? "a" : "r"

missingval(var::CDM.AbstractVariable, md::Metadata{<:CDMsource}) =
    missingval(md)
missingval(var::CDM.AbstractVariable, args...) = 
    missingval(Metadata{sourcetrait(var)}(CDM.attribs(var)))

@inline function get_scale(metadata::Metadata{<:CDMsource}, scaled::Bool)
    scale = scaled ? get(metadata, "scale_factor", nothing) : nothing
    offset = scaled ? get(metadata, "add_offset", nothing) : nothing
    return scale, offset
end

# Rasters methods for CDM types ###############################

Raster(ds::AbstractVariable; kw...) = _raster(ds; kw...)

function _open(f, source::CDMsource, filename::AbstractString; write=false, kw...)
    checkfilename(source, filename)
    ds = sourceconstructor(source)(filename, openmode(write))
    _open(f, source, ds; kw...)
end
function _open(f, source::CDMsource, ds::AbstractDataset; 
    name=nokw, 
    group=nothing, 
    mod=NoMod(), 
    kw...
)
    g = _getgroup(ds, group)
    if isnokw(name)
        cleanreturn(f(g)) 
    else
        key = string(_name_or_firstname(g, name))
        v = CDM.variable(g, key)
        _open(f, source, v; mod)
    end
end
_open(f, ::CDMsource, var::AbstractArray; mod=NoMod(), kw...) = 
    cleanreturn(f(_maybe_modify(var, mod)))

# This allows arbitrary group nesting
_getgroup(ds, ::Union{Nothing,NoKW}) = ds
_getgroup(ds, group::Union{Symbol,AbstractString}) = ds.group[String(group)]
_getgroup(ds, group::Pair) = _getgroup(ds.group[String(group[1])], group[2])

filename(ds::AbstractDataset) = CDM.path(ds)
filename(var::AbstractVariable) = CDM.path(CDM.dataset(var))

filekey(ds::AbstractDataset, name::Union{String,Symbol}) = Symbol(name)
filekey(ds::AbstractDataset, name) = _name_or_firstname(ds, name)
cleanreturn(A::AbstractVariable) = Array(A)
haslayers(::CDMsource) = true
defaultcrs(::CDMsource) = EPSG(4326)
defaultmappedcrs(::CDMsource) = EPSG(4326)

function _dims(var::AbstractVariable{<:Any,N}, crs=nokw, mappedcrs=nokw) where N
    dimnames = CDM.dimnames(var)
    ntuple(Val(N)) do i
        _cdm_dim(CDM.dataset(var), dimnames[i], crs, mappedcrs)
    end
end
_metadata(var::AbstractVariable; attr=CDM.attribs(var)) =
    _metadatadict(sourcetrait(var), attr)

function _dims(ds::AbstractDataset, dimdict::Dict)
    map(CDM.dimnames(ds)) do dimname
        dimdict[dimname]
    end |> Tuple
end
_metadata(ds::AbstractDataset; attr=CDM.attribs(ds)) =
    _metadatadict(sourcetrait(ds), attr)
function _layerdims(ds::AbstractDataset; layers, dimdict)
    map(layers.vars) do var
        map(CDM.dimnames(var)) do dimname
            basedims(dimdict[dimname])
        end |> Tuple
    end
end
function _layermetadata(ds::AbstractDataset; layers)
    map(layers.attrs) do attr
        md = _metadatadict(sourcetrait(ds), attr)
        if haskey(attr, "grid_mapping")
            if haskey(ds, attr["grid_mapping"])
                md["grid_mapping"] = Dict(CDM.attribs(ds[attr["grid_mapping"]]))
            else
                global_attrs = CDM.attribs(ds)
                if haskey(global_attrs, attr["grid_mapping"])
                    md["grid_mapping"] = global_attrs["grid_mapping"]
                end
            end     
        end
        md
    end
end


# Utils ########################################################################

# Get the names of all data layers, removing
# coordinate, bounds, and geometry variables etc
function _layer_names(ds) 
    names = keys(ds)
    return _layer_names(ds, names, Dict{String,Any}(), Dict{String,Any}())
end
function _layer_names(ds, var_names, vars_dict, attrs_dict)
    bounds_names = String[]
    geometry_names = String[]
    grid_mapping_names = String[]
    coordinate_names = String[]
    node_coordinate_names = String[]
    # Sometimes coordinates are not included in var attributes
    # So we check standard names to know vars are lookups and not data
    standard_names = String[]
    for n in zip(var_names)
        get!(() -> CDM.variable(ds, n), vars_dict, n)
        get!(() -> CDM.attribs(vars_dict[n]), attrs_dict, n)
        # Coordinate variables
        haskey(attr, "coordinates") && setdiff(coordinate_names, _cdm_coordinates(attr))
        haskey(attr, "node_coordinates") && setdiff(node_coordinate_names, _cdm_node_coordinates(attr))
        # Bounds variables
        haskey(attr, "bounds") && push!(bounds_names, attr["bounds"])
        # Geometry variables
        haskey(attr, "geometry") && push!(geometry_names, attr["geometry"])
        # Grid mapping variables
        haskey(attr, "grid_mapping") && push!(grid_mapping_names, attr["grid_mapping"])
        # Common Standard names - very likely coordinates
        haskey(attr, "standard_name") && attr["standard_name"] in keys(CDM_STANDARD_NAME_MAP) && 
            push!(standard_names, n)
    end
    dim_names = CDM.dimnames(ds)
    # if all(n -> n in var_names, dim_names)
    # end
    return setdiff(varnames, bounds_names, geometry_names, grid_mapping_names, coordinate_names, node_coordinate_names, standard_names)
end

_cdm_coordinates(attr) = collect(split(attr["coordinates"], ' '))
_cdm_node_coordinates(attr) = collect(split(attr["node_coordinates"], ' '))

function _organise_dataset(ds::AbstractDataset, names=nokw, group::NoKW=nokw)
    var_names = keys(ds)
    vars_dict = Dict{String,Any}()
    attrs_dict = Dict{String,Any}()
    layer_names = isnokw(names) ? _layer_names(ds, var_names, vars_dict, attrs_dict) : names

    dim_dict = Dict{String,Dimension}()
    geometry_dict = Dict{String,Any}()
    grid_mapping_dict = Dict{String,Any}()
    output_layers_vec = Vector{CDM.AbstractVariable}(undef, length(layer_names))
    output_attrs_vec = Vector{CDM.Attributes}(undef, length(layer_names))

    # Loop over the layers we want to load as rasters
    # As we go, we add dimensions to dim_dict to be used in other layers
    for (i, var_name) in enumerate(layer_names)
        var = get(() -> CDM.variable(ds, var_name), vars_dict, var_name)
        attr = get(() -> CDM.attribs(var), attrs_dict, var_name)
        coord_names = _cdm_coordinates(attr)
        layer_dimnames_vec[i] = dimnames = CDM.dimnames(var)
        attr = attrs_dict[var_name]
        geometry_key = if haskey(attr, "geometry")
            attr["geometry"]
        end
        grid_mapping_key = if haskey(attr, "grid_mapping")
            attr["grid_mapping"]
        end
        # Loop over the dimensions of this layer, adding missing dims to dim_dict 
        for dimname in dimnames 
            haskey(dim_dict, dimname) && continue
            dim_dict[dimname] = _cdm_dim(
                ds, vars_dict, attrs_dict, geometry_dict, grid_mapping_dict, geometry_key, grid_mapping_key, coord_names, dimname
            )
        end
        output_layers_vec[i] = var
        output_attrs_vec[i] = attr
    end
    # TODO output
    (; layers, dims, layerdims_vec, layermetadata_vec)
end
_organise_dataset(ds::AbstractDataset, names, group) =
    _organise_dataset(ds.group[group], names, nokw)

function _cdm_dim(ds, vars_dict, attrs_dict, geometry_dict, grid_mapping_dict, geometry_key, grid_mapping_key, coord_names, dimname)
    crs = if !isnothing(grid_mapping_key) 
        get(grid_mapping_dict, grid_mapping_key) do
            _cdm_grid_mapping(ds, grid_mapping_key)
        end
    end
    if haskey(ds, dimname)
        # The dimension has a matching variable
        dim_attrs = attrs_dict[dimname]
        D = _cdm_dimtype(dim_attrs, dimname)
        lookup = _cdm_lookup(ds, dimname, D, dim_attrs)
        D(lookup)
    else
        if !isnothing(geometry_key)
            geom = get(geometry_dict, geometry_key) do
                _cdm_geometry(ds, geometry_key)
            end
            # Check that geometry node coordinates are for this dimension
            is_geometry_dimension = all(_cdm_node_coordinates(geom)) do nc
                dimname in CDM.dimnames(vars_dict[nc])
            end
            if is_geometry_dimension
                D = _cdm_dimtype(NoMetadata(), d)
                lookup = _cdm_geometry_lookup(vars_dicts, attrs_dict, geom, dimname)
                return D(lookup)
            end
        end
        # No matching variable. Get all the coordinates for this dimension
        dimension_coord_names = filter(coordnames) do coordname
            dimnames = CDM.dimnames(vars_dict[coordname])
            length(dimnames) == 1 && only(dimnames) == dimname  
        end
        # If there are coordinates for this dimension
        if !isempty(dimension_coord_names)
            # TODO rotations
            dims = map(dimension_coord_names) do d
                dim_attrs = attrs_dict[d]
                D = _cdm_dimtype(dim_attrs, d)
                lookup = _cdm_lookup(ds, d, D, dim_attrs)
                D(lookup)
            end
            D = _cdm_dimtype(NoMetadata(), dimname)
            # Use a MergedLookup with all dimensions
            mergedims(dims, dims => D)
        else
            # The var doesn't exist. Maybe its `complex` or some other marker,
            # so make it a custom `Dim` with `NoLookup`
            lookup = NoLookup(Base.OneTo(len))
            D = _cdm_dimtype(NoMetadata(), dimname)
            D(lookup)
        end
    end
end

# TODO don't load all keys here with _layers
_name_or_firstname(ds::AbstractDataset, name) = Symbol(name)
function _name_or_firstname(ds::AbstractDataset, name::Union{Nothing,NoKW}=_layer_names(ds))
    if length(names) > 0
        return Symbol(first(names))
    else
        throw(ArgumentError("No non-dimension layers found in dataset with keys: $(keys(ds))"))
    end
end

# Find the matching dimension constructor. If its an unknown name
# use the generic Dim with the dim name as type parameter
function _cdm_dimtype(attrib, dimname)
    # First use the axis - X/Y/Z/T as the 
    # closest match to DimensionalData standard dimension names
    if haskey(attrib, "axis")
        k = attrib["axis"]
        if haskey(CDM_AXIS_MAP, k)
            return CDM_AXIS_MAP[k]
        end
    end
    # Then use standard_name - latitude/longitude/height/time etc
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

# _cdm_lookup
# Generate a `Lookup` from a CDM dim.
function _cdm_lookup(ds::AbstractDataset, attr, dimname, D::Type, attr)
    var = ds[dimname]
    data = Missings.disallowmissing(var[:])
    metadata = _metadatadict(sourcetrait(ds), attr)
    return _cdm_lookup(data, ds, var, attr, dimname, D, metadata)
end
# For unknown types we just make a Categorical lookup
function _cdm_lookup(data::AbstractArray, ds::AbstractDataset, var, attr, dimname, D::Type, metadata)
    Categorical(data; order=Unordered(), metadata=metadata)
end
# For Number and AbstractTime we generate order/span/sampling
# We need to include `Missing` in unions in case `_FillValue` is used
# on coordinate variables in a file and propagates here.
function _cdm_lookup(
    data::AbstractArray{<:Union{Missing,Number,Dates.AbstractTime}},
    ds::AbstractDataset, var, attr, dimname, D::Type, metadata,
)
    # Assume the locus is at the center of the cell if boundaries aren't provided.
    # http://cfconventions.org/cf-conventions/cf-conventions.html#cell-boundaries
    order = LA.orderof(index)
    # Detect lat/lon
    span, sampling = if eltype(index) <: Union{Missing,Dates.AbstractTime}
        _cdm_period(index, metadata)
    else
        _cdm_span(index, order)
    end
    # We only use Explicit if the span is not Regular
    # This is important for things like rasterizatin and conversion 
    # to gdal to be easy, and selectors are faster.
    # TODO are there any possible floating point errors from this?
    if haskey(attr, "bounds")
        span, sampling = if isregular(span)
            span, Intervals(Center())
        else
            boundskey = attr["bounds"]
            boundsmatrix = Array(ds[boundskey])
            locus = if mapreduce(==, &, view(boundsmatrix, 1, :), index)
                Start()
            elseif mapreduce(==, &, view(boundsmatrix, 2, :), index)
                End()
            else
                # TODO better checks here
                Center()
            end
            Explicit(boundsmatrix), Intervals(locus)
        end
    end
    return _cdm_lookup(data, D, order, span, sampling, metadata)
end
# For X and Y use a Mapped <: AbstractSampled lookup
function _cdm_lookup(
    data, D::Type{<:Union{<:XDim,<:YDim}}, order::Order, span, sampling, metadata
)
    # If units are degrees north/east, use EPSG:4326 as the mapped crs
    units = get(metadata, "units", "")
    mappedcrs = if isnokw(mappedcrs) && (units in CF_DEGREES_NORTH || units in CF_DEGREES_EAST)
        EPSG(4326)
    end
    # Additionally, crs and mappedcrs should be identical for Regular lookups
    crs = if span isa Regular
        mappedcrs
    end
    dim = DD.basetypeof(D)()
    return Mapped(data; order, span, sampling, metadata, crs, mappedcrs, dim)
end
# Band dims have a Categorical lookup, with order
# This is not CF, just for consistency with GDAL
function _cdm_lookup(data, D::Type{<:Band}, order::Order, span, sampling, metadata)
    Categorical(data, order, metadata)
end
# Otherwise use a Sampled lookup
function _cdm_lookup(data, D::Type, order::Order, span, sampling, metadata)
    Sampled(data, order, span, sampling, metadata)
end

function _cdm_span(index, order)
    # Handle a length 1 index
    length(index) == 1 && return Regular(zero(eltype(index))), Points()

    step = if eltype(index) <: AbstractFloat
        # Calculate step, avoiding as many floating point errors as possible
        st = Base.step(Base.range(Float64(first(index)), Float64(last(index)); length = length(index)))
        st_rd = round(st, digits = Base.floor(Int,-log10(eps(eltype(index))))) # round to nearest digit within machine epsilon
        isapprox(st_rd, st; atol = eps(eltype(index))) ? st_rd : st # keep the rounded number if it is very close to the original
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
function _cdm_period(index, metadata::Metadata{<:CDMsource})
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

function _attribdict(md::Metadata{<:CDMsource}) 
    attrib = Dict{String,Any}()
    for (k, v) in md
        v isa Tuple && continue
        attrib[string(k)] = v
    end
    return attrib
end
_attribdict(md) = Dict{String,Any}()

# Add axis and standard name attributes to dimension variables
# We need to get better at guaranteeing if X/Y is actually measured in `longitude/latitude`
# CF standards requires that we specify "units" if we use these standard names
_cdm_set_axis_attrib!(atr, dim::X) = atr["axis"] = "X" # at["standard_name"] = "longitude";
_cdm_set_axis_attrib!(atr, dim::Y) = atr["axis"] = "Y" # at["standard_name"] = "latitude";
_cdm_set_axis_attrib!(atr, dim::Z) = (atr["axis"] = "Z"; atr["standard_name"] = "depth")
_cdm_set_axis_attrib!(atr, dim::Ti) = (atr["axis"] = "T"; atr["standard_name"] = "time")
_cdm_set_axis_attrib!(atr, dim) = nothing

_cdm_shiftlocus(dim::Dimension) = _cdmshiftlocus(lookup(dim), dim)
_cdm_shiftlocus(::Lookup, dim::Dimension) = dim
function _cdm_shiftlocus(lookup::AbstractSampled, dim::Dimension)
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

function Base.write(filename::AbstractString, source::CDMsource, A::AbstractRaster;
    append=false,
    force=false,
    kw...
)
    mode = checkwritemode(source, filename, append, force)
    ds = sourceconstructor(source)(filename, mode; attrib=_attribdict(metadata(A)))
    try
        writevar!(ds, source, A; kw...)
    finally
        close(ds)
    end
    return filename
end
function Base.write(filename::AbstractString, source::Source, s::AbstractRasterStack{K,T};
    append=false,
    force=false,
    missingval=nokw,
    f=identity,
    kw...
) where {Source<:CDMsource,K,T}
    mode = checkwritemode(source, filename, append, force)
    ds = sourceconstructor(source)(filename, mode; attrib=_attribdict(metadata(s)))
    missingval = _stack_nt(s, isnokw(missingval) ? Rasters.missingval(s) : missingval)
    try
        mods = map(keys(s)) do k
            writevar!(ds, source, s[k]; missingval=missingval[k], kw...)
        end
        f(OpenStack{Source,K,T}(ds, mods))
    finally
        close(ds)
    end
    return filename
end

# Add a var array to a dataset before writing it.
function writevar!(ds::AbstractDataset, source::CDMsource, A::AbstractRaster{T,N};
    verbose=true,
    missingval=nokw,
    metadata=nokw,
    chunks=nokw,
    chunksizes=_chunks_to_tuple(A, dims(A), chunks),
    scale=nokw,
    offset=nokw,
    coerce=convert,
    eltype=Missings.nonmissingtype(T),
    write=true,
    name=DD.name(A),
    options=nokw,
    driver=nokw,
    f=identity,
    kw...
) where {T,N}
    _check_allowed_type(source, eltype)
    write = f === identity ? write : true
    _def_dim_var!(ds, A)
    metadata = if isnokw(metadata) 
        DD.metadata(A)
    elseif isnothing(metadata)
        NoMetadata()
    else
        metadata
    end

    missingval_pair = _write_missingval_pair(A, missingval; eltype, verbose, metadata)

    attrib = _attribdict(metadata)
    # Scale and offset
    scale = if isnokw(scale) || isnothing(scale)
        delete!(attrib, "scale_factor")
        nothing
    else
        attrib["scale_factor"] = scale
    end
    offset = if isnokw(offset) || isnothing(offset)
        delete!(attrib, "add_offset")
        nothing
    else
        attrib["add_offset"] = offset
    end

    mod = _mod(eltype, missingval_pair, scale, offset, coerce)

    if !isnothing(missingval_pair[1])
        attrib["_FillValue"] = missingval_pair[1]
    end

    key = if isnokw(name) || string(name) == ""
        UNNAMED_FILE_KEY
    else
        string(name)
    end

    dimnames = lowercase.(string.(map(Rasters.name, dims(A))))
    var = CDM.defVar(ds, key, eltype, dimnames; attrib, chunksizes, kw...)

    if write
        m = _maybe_modify(parent(var), mod)
        # Write with a DiskArays.jl broadcast
        m .= A
        # Apply `f` while the variable is open
        f(m)
    end

    return mod
end

_check_allowed_type(trait, eltyp) = nothing
function _check_allowed_type(::CDMsource, eltyp)
    eltyp <: CDMallowedType || throw(ArgumentError("""
    Element type $eltyp cannot be written to NetCDF. Convert it to one of $(Base.uniontypes(CDMallowedType)),
    usually by broadcasting the desired type constructor over the `Raster`, e.g. `newrast = Float32.(rast)`"))
    """
    ))
end

_def_dim_var!(ds::AbstractDataset, A) = map(d -> _def_dim_var!(ds, d), dims(A))
function _def_dim_var!(ds::AbstractDataset, dim::Dimension)
    dimname = lowercase(string(DD.name(dim)))
    haskey(ds.dim, dimname) && return nothing
    CDM.defDim(ds, dimname, length(dim))
    _def_lookup_var!(ds, dim, dimname)
end 
_def_lookup_var!(ds::AbstractDataset, dim::Dimension{<:AbstractNoLookup}, dimname) = nothing
function _def_lookup_var!(ds::AbstractDataset, dim::Dimension{<:ArrayLookup}, dimname)  
    # TODO what do we call the geometry variable, what if more than 1
    return nothing
end
function _def_lookup_var!(ds::AbstractDataset, dim::Dimension{<:GeometryLookup}, dimname)  
    l = lookup(dim)
    geom = _geometry_cf_encode(parent(l))
    # Define base geometry attributes

    coord_names_raw = map(dims(l)) do d
        string(DD.name(d))
    end

    coord_names = if any(n -> haskey(ds, n), coord_names_raw)
        (dimname,) .* coord_names_raw
    else
        coord_names_raw
    end

    attrib = Dict{String,Any}(
        "geometry_type" => geom.geometry_type,
        "node_coordinates" => join(coord_names, " "),
    )
    if haskey(geom, :node_count)
        # More nodes than dimension length: define a node dimension
        node_dimname = dimname * "_node"
        CDM.defVar(ds, "node_count", geom.node_count, dimname)
        CDM.defDim(ds, node_dimname, length(geom.x))
        attrib["node_count"] = dimname * "_node_count"
    else
        # Node count is dimension length
        node_dimame = dimname
    end
    if haskey(geom, :part_node_count)
        # We also have geometry parts
        part_dimname = dimname * "part"
        part_varname = dimname * "_part_node_count"
        CDM.defDim(ds, part_dimname, length(geom.part_node_count))
        CDM.defVar(ds, part_varname, geom.part_node_count, part_dimname)
        attrib["part_node_count"] = part_varname
    end
    if haskey(geom, :interior_ring)
        # And interior rings
        ring_varname = dimname * "_interior_ring"
        CDM.defVar(ds, ring_varname, geom.interior_ring, part_dimname)
        attrib["interior_ring"] = ring_varname
    end
    # Define coordinate variables
    CDM.defVar(ds, coord_names[1], geom.x, node_dimame)
    CDM.defVar(ds, coord_names[2], geom.y, node_dimame)
    # TODO: add z and m, if present
    # And the geometry container
    CDM.defVar(ds, "geometry_container", nothing, (), attrib)
    return nothing
end
function _def_lookup_var!(ds::AbstractDataset, dim::Dimension{<:Union{Sampled,AbstractCategorical}}, dimname)  
    # Shift lookup before conversion to Mapped
    dim = _cdm_shiftlocus(dim)
    if dim isa Y || dim isa X
        dim = convertlookup(Mapped, dim)
    end
    # Attributes
    attrib = _attribdict(metadata(dim))
    _cdm_set_axis_attrib!(attrib, dim)
    # Bounds variables
    if sampling(dim) isa Intervals
        bounds = Dimensions.dim2boundsmatrix(dim)
        boundskey = get(metadata(dim), :bounds, string(dimname, "_bnds"))
        push!(attrib, "bounds" => boundskey)
        CDM.defVar(ds, boundskey, bounds, ("bnds", dimname))
    end
    CDM.defVar(ds, dimname, Vector(index(dim)), (dimname,); attrib)
    return nothing
end

function _classify_variables(ds::AbstractDataset)
end


function _read_geometry(ds, key)
    var = CDM.variable(ds, key)
    var_attrib = CDM.attribs(var)
    geom_key = var_attrib["geometry"]
    geom_var = CDM.variable(ds, geom_key)
    attrib = CDM.attribs(geom_var)
    geometry_type = attrib["geometry_type"]

    node_coordinate_names = split(attrib["node_coordinates"], ' ')
    node_coordinate_vars = map(c -> CDM.variable(ds, c), node_coordinate_names)
    # Map the node coordinates to a vector of n-tuples.
    node_coordinates = map(tuple, node_coordinate_vars...)
    # There are three types of geometries:
    # - point (can be point or multipoint depending on presence of node_count)
    # - line (can be linestring or multilinestring depending on presence of part_node_count)
    # - polygon (can be poly or multipoly depending on presence of part_node_count, multipoly also needs interior_ring)
    # This is the meat of the function, that parses the geometry and returns a named tuple 
    trait, geom_parsing_args = if geometry_type == "point"
        if haskey(attrib, "node_count") # actually multipoint
            (GI.MultiPointTrait(), (; node_coordinates, node_count = CDM.variable(ds, attrib["node_count"])))
        else # regular points
            (GI.PointTrait(), (; node_coordinates))
        end
    elseif geometry_type == "line"
        if haskey(attrib, "part_node_count") # actually multilinestring
            (GI.MultiLineStringTrait(), (; 
                node_coordinates, 
                node_count = CDM.variable(ds, attrib["node_count"]), 
                part_node_count = CDM.variable(ds, attrib["part_node_count"])
            ))
        else # regular linestring
            (GI.LineStringTrait(), (; 
                node_coordinates, 
                node_count = CDM.variable(ds, attrib["node_count"])
            ))
        end
    elseif geometry_type == "polygon"
        if haskey(attrib, "part_node_count") # actually multipoly
            if haskey(attrib, "interior_ring") # this will usually be the case
                (GI.MultiPolygonTrait(), (; 
                    node_coordinates, 
                    node_count = CDM.variable(ds, attrib["node_count"]), 
                    part_node_count = CDM.variable(ds, attrib["part_node_count"]), 
                    interior_ring = CDM.variable(ds, attrib["interior_ring"])
                ))
            else # no interior ring, but part node count is present
                # so, we need to assert that part_node_count has the same length as node_count
                # in which case, interior_ring is always false by definition.
                node_count = CDM.variable(ds, attrib["node_count"])
                part_node_count = CDM.variable(ds, attrib["part_node_count"])

                if length(node_count) != length(part_node_count)
                    throw(ArgumentError("""
                    When parsing geometry lookup for CDM variable $key, 
                    the length of `node_count` ($(length(node_count))) and 
                    `part_node_count` ($(length(part_node_count))) must be the same.
                    """))
                end
                (GI.PolygonTrait(), (; node_coordinates, node_count))
            end
        else # regular polygon
            (GI.PolygonTrait(), (; node_coordinates, node_count = CDM.variable(ds, attrib["node_count"])))
        end
    else
        throw(ArgumentError("""
        When parsing geometry lookup for CDM variable $key, 
        the geometry type is $geometry_type,
        but we could not figure out what that is!  
        Allowed geometry types are: `point`, `line`, `polygon`.

        If you need something else, please open a Github issue!
        """))
    end
    
    crs = if haskey(attrib, "grid_mapping")
        grid_mapping_var = CDM.variable(ds, attrib["grid_mapping"])
        grid_mapping_attrib = CDM.attribs(grid_mapping_var)
        if haskey(grid_mapping_attrib, "crs_wkt")
            WellKnownText2(GeoFormatTypes.CRS(), grid_mapping_attrib["crs_wkt"])
        elseif haskey(grid_mapping_attrib, "spatial_epsg")
            EPSG(parse(Int, grid_mapping_attrib["spatial_epsg"]))
        elseif haskey(grid_mapping_attrib, "proj4string")
            ProjString(grid_mapping_attrib["proj4string"])
        else
            nothing # unparseable - TODO get CFProjections.jl to do it.
        end
    else
        nothing
    end

    # return the geometry directly
    # TODO: parse coordinate names, to get dimensions
    # depends on what utils / luts we have elsewhere.
    return _geometry_cf_decode(trait, ds, geom_parsing_args; crs)
    
    # (; dimension, internal_dimensions, geometry_type, node_count, part_node_count, interior_ring, node_coordinates, grid_mapping)
end