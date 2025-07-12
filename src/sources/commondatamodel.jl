const CDM = CommonDataModel

const CDMallowedType = Union{Int8,UInt8,Int16,UInt16,Int32,UInt32,Int64,UInt64,Float32,Float64,Char,String}

const CF_DEGREES_NORTH = ("degrees_north", "degree_north", "degree_N", "degrees_N", "degreeN", "degreesN")
const CF_DEGREES_EAST = ("degrees_east", "degree_east", "degree_E", "degrees_E", "degreeE", "degreesE")

const UNNAMED_FILE_KEY = "unnamed"

const CDM_PERIODS = Dict(
    "years" => Year,
    "months" => Month,
    "days" => Day,
    "hours" => Hour,
    "minutes" => Minute,
    "seconds" => Second,
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

_cf_name(dim::Dimension) = lowercase(string(name(dim)))
_cf_name(dim::Ti) = "time"
_cf_name(dim::X) = "x"
_cf_name(dim::Y) = "y"
_cf_name(dim::Z) = "z"

_cf_axis(dim::Dimension) = nothing
_cf_axis(dim::Ti) = "T"
_cf_axis(dim::X) = "X"
_cf_axis(dim::Y) = "Y"
_cf_axis(dim::Z) = "Z"

_cf_shortname(dim::Dimension) = nothing
_cf_shortname(dim::Ti) = "time"
_cf_shortname(dim::X) = "lon"
_cf_shortname(dim::Y) = "lat"
_cf_shortname(dim::Z) = "height"

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

missingval(var::AbstractArray, md::Metadata{<:CDMsource}) = missingval(md)
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

# Converts Char arrays to String folowing CF chapter 2
struct DiskCharToString{N,M} <: DiskArrays.AbstractDiskArray{String,N}
    parent::AbstractArray{Char,M} 
end
DiskCharToString(A::AbstractArray{Char,M}) where M = DiskCharToString{M-1,M}(A)

Base.parent(A::DiskCharToString) = A.parent
Base.size(A::DiskCharToString) = size(parent(A))[2:end]
DiskArrays.haschunks(A::DiskCharToString) = DiskArrays.haschunks(parent(A))
function DiskArrays.readblock!(A::DiskCharToString, dest, I...)
    src = parent(A)[:, I...]
    for I in CartesianIndices(dest)
        dest[I] = _chars_to_string(view(src, :, I))
    end
    return dest
end

_chars_to_string(chars) = String(strip(join(chars), '\0'))

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

function _dims(var::AbstractVariable{<:Any,N}, crs=nokw, mappedcrs=nokw) where N
    dimnames = CDM.dimnames(var)
    ntuple(Val(N)) do i
        _cdm_dim(CDM.dataset(var), dimnames[i], crs, mappedcrs)
    end
end
_metadata(var::AbstractVariable; attr=CDM.attribs(var)) =
    _metadatadict(sourcetrait(var), attr)

_metadata(ds::AbstractDataset; attr=CDM.attribs(ds)) =
    _metadatadict(sourcetrait(ds), attr)
function _layerdims(ds::AbstractDataset; layers, dimdict)
    map(layers.vars) do var
        map(CDM.dimnames(var)) do dimname
            basedims(dimdict[dimname])
        end |> Tuple
    end
end
function _layermetadata(ds::AbstractDataset; attrs)
    map(attrs) do attr
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


# Get the names of all data layers, removing
# coordinate, bounds, and geometry variables etc
function _layer_names(ds) 
    names = keys(ds)
    return _layer_names(ds, names, Dict{String,Any}(), Dict{String,Any}())
end
function _layer_names(ds, var_names, vars_dict, attrs_dict)
    # Keep everything separate for clarity
    bounds_names = String[]
    geometry_names = String[]
    grid_mapping_names = String[]
    coordinate_names = String[]
    units_metadata_names = String[]
    node_coordinate_names = String[]
    quantization_names = String[]
    climatology_names = String[]
    # Sometimes coordinates are not included in var attributes
    # So we check standard names to know vars are lookups and not data
    standard_names = String[]
    for n in var_names
        var = _get_var!(vars_dict, ds, n)
        attr = _get_attr!(attrs_dict, vars_dict, ds, n)
        # Coordinate variables
        haskey(attr, "coordinates") && append!(coordinate_names, collect(_cdm_coordinates(attr)))
        haskey(attr, "node_coordinates") && append!(node_coordinate_names, collect(_cdm_node_coordinates(attr)))
        # Bounds variables
        haskey(attr, "bounds") && push!(bounds_names, attr["bounds"])
        # Quantization variables (8.4.2)
        haskey(attr, "quantization") && push!(quantization_names, attr["quantization"])
        # Climatology variables (7.4)
        haskey(attr, "climatology") && push!(climatology_names, attr["climatology"])
        # Geometry variables (7.5)
        if haskey(attr, "geometry") 
            geometry_container_name = attr["geometry"]
            push!(geometry_names, geometry_container_name)
            geometry_attr = _get_attr!(attrs_dict, vars_dict, ds, geometry_container_name)
            haskey(geometry_attr, "node_coordinates") && push!(geometry_names, _cf_node_coordinate_names(geometry_attr)...)
            haskey(geometry_attr, "node_count") && push!(geometry_names, geometry_attr["node_count"])
            haskey(geometry_attr, "part_node_count") && push!(geometry_names, geometry_attr["part_node_count"])
            haskey(geometry_attr, "interior_ring") && push!(geometry_names, geometry_attr["interior_ring"])
        end
        # Grid mapping variables - there may be multiple
        haskey(attr, "grid_mapping") && append!(grid_mapping_names, _grid_mapping_keys(attr["grid_mapping"]))
        # Units metadata variables
        haskey(attr, "units_metadata") && push!(units_metadata_names, attr["units_metadata"])
        # Common Standard names - very likely coordinates, but is this needed?
        # haskey(attr, "standard_name") && attr["standard_name"] in keys(CDM_STANDARD_NAME_MAP) && 
            # push!(standard_names, n)
    end
    dim_names = CDM.dimnames(ds)
    return setdiff(var_names,
        dim_names,
        bounds_names,
        quantization_names,
        climatology_names,
        geometry_names,
        grid_mapping_names,
        coordinate_names,
        units_metadata_names,
        node_coordinate_names,
        standard_names
    )
end

function _cdm_coordinates(attr)
    if haskey(attr, "coordinates")
        String.(split(attr["coordinates"], ' '))
    else
        String[]
    end
end
function _cdm_node_coordinates(attr)
    if haskey(attr, "node_coordinates")
        String.(split(attr["node_coordinates"], ' '))
    else
        String[]
    end
end

#= Here we convert the surface of the CF standard that we support 
into DimensionalData.jl objects
Unhandled
- ancillary variables 3.4 (connection ignored)
- flag values 3.5 (ignored, raw data used - could use CategoricalArrays here?)
- parametric vertical coordinates 4.3.3 (ignored)
=# 
function _organise_dataset(ds::AbstractDataset, names=nokw, group::NoKW=nokw)
    # Start with the names of all variables
    var_names = keys(ds)
    # Define vars and attrs dicts so these are only loaded once
    vars_dict = Dict{String,Any}()
    attrs_dict = Dict{String,Any}()
    geometry_dict = Dict{String,Any}()
    crs_dict = Dict{String,Any}()
    # Define a dimensions dict so these are only generated once
    dim_dict = Dict{String,Dimension}()
    # Refdims need a consistent order as it is not changed later
    refdim_dict = OrderedDict{String,Dimension}()
    # Get the layer names to target
    layer_names = isnokw(names) ? _layer_names(ds, var_names, vars_dict, attrs_dict) : names
    layer_dimnames_vec = Vector{Tuple}(undef, length(layer_names))

    # Define output data and metadata vars
    output_layers_vec = Vector{AbstractArray}(undef, length(layer_names))
    output_layerdims_vec = Vector{Tuple}(undef, length(layer_names))
    output_attrs_vec = Vector{Any}(undef, length(layer_names))

    # Loop over the layers we want to load as rasters
    # As we go, we add dimensions to dim_dict to be used in other layers
    for (i, var_name) in enumerate(layer_names)
        # Get the variable and its attributes
        var = get(() -> CDM.variable(ds, var_name), vars_dict, var_name)
        attr = get(() -> CDM.attribs(var), attrs_dict, var_name)
        # Get its dimensions
        layer_dimnames_vec[i] = dimnames = CDM.dimnames(var)
        # Get its coordinates
        coord_names = _cdm_coordinates(attr)
        # Remove coordinates from metadata
        if haskey(attr, "coordinates")
            delete!(attr, "coordinates")
        end
        # Check for a geometry variable (CF 7.5)
        geometry_key = if haskey(attr, "geometry")
            pop!(attr, "geometry")
        end
        # Check for a grid mapping variable
        if haskey(attr, "grid_mapping")
            grid_mapping_key = pop!(attr, "grid_mapping")
            crs = get!(() -> _cf_crs(ds, grid_mapping_key), crs_dict, grid_mapping_key)
        else
            crs = grid_mapping_key = nothing
        end
        if eltype(var) <: Char && ndims(var) > 1 
            # Remove the char dimension
            layer_dimnames_vec[i] = dimnames = dimnames[2:end]
        end
        # Loop over the dimensions of this layer, adding missing dims to dim_dict 
        output_layerdims_vec[i] = map(dimnames) do dimname
            get!(dim_dict, dimname) do
                _cdm_dim(
                    ds, vars_dict, attrs_dict, geometry_dict, geometry_key, var_name, coord_names, dimname, crs, grid_mapping_key
                )
            end
        end
        unformatted_dims = map(dimnames) do dimname
            dim_dict[dimname]
        end
        formatted_dims = format(unformatted_dims, var)
        map(dimnames, formatted_dims) do dimname, d
            dim_dict[dimname] = d
        end
        # Find scalar coordinates to use as refdims (5.7)
        for coord_name in coord_names
            haskey(refdim_dict, coord_name) && continue
            coord_var = _get_var!(vars_dict, ds, coord_name)
            if ndims(coord_var) == 0
                coord_attr = _get_attr!(attrs_dict, vars_dict, ds, coord_name)
                data = reshape(collect(coord_var), 1)
                D = _cdm_dimtype(coord_attr, coord_name)
                metadata = _metadatadict(sourcetrait(ds), attr)
                lookup = _cdm_lookup(data, ds, coord_var, coord_attr, coord_name, D, metadata, nothing)
                refdim = D(lookup)
                refdim_dict[coord_name] = refdim
            end
        end
        # Add this layer to the output
        output_layers_vec[i] = var
        output_attrs_vec[i] = attr
    end

    return (; 
        names_vec=layer_names, 
        layers_vec=output_layers_vec, 
        layerdims_vec=output_layerdims_vec, 
        layermetadata_vec=output_attrs_vec, 
        dim_dict,
        refdim_dict,
    )
end
_organise_dataset(ds::AbstractDataset, names, group) =
    _organise_dataset(ds.group[group], names, nokw)


function _cdm_dim(ds, vars_dict, attrs_dict, geometry_dict, geometry_key, var_name, coord_names, dimname, crs, grid_mapping_key)
    # First check for geometry dimensions (remove most complicted first)
    if is_geometry(ds, vars_dict, attrs_dict, geometry_key, dimname)
        return _geometry_dim(ds, vars_dict, attrs_dict, geometry_dict, geometry_key, crs)
    end

    # Then get all the coordinates for this dimension
    dimension_coord_names = filter(coord_names) do coordname
        coord_dimnames = CDM.dimnames(_get_var!(vars_dict, ds, coordname))
        dimname in coord_dimnames
    end
    # Now dimension/lookup will depend on if there are associated coordinates 
    if isempty(dimension_coord_names)
        if haskey(ds, dimname) # The var name is the dim name
            return _standard_dim(ds, vars_dict, attrs_dict, dimname, crs)
        else # A matching variable doesn't exist. 
            return _discrete_axis_dim(ds, dimname)
        end
    else # Use coord names
        # Otherwise dims/lookup depend on the number of associated coords
        if length(dimension_coord_names) == 1
            if haskey(ds, dimname) # The var name is the dim name
                # Only one coordinate, so just use it as the dimension/lookup
                return _standard_dim(ds, vars_dict, attrs_dict, dimname, crs)
            else
                return _standard_dim(ds, vars_dict, attrs_dict, dimension_coord_names[1], crs)
            end
        else
            first_coord_name = dimension_coord_names[1]
            first_coord_var = _get_var!(vars_dict, ds, first_coord_name)
            # Short circuit for char dims
            if is_char_categorical(first_coord_var)
                first_coord_attrs = _get_attr!(attrs_dict, vars_dict, ds, first_coord_name)
                return _char_categorical_dim(ds, first_coord_var, first_coord_attrs, first_coord_name)
            end
            if is_multidimensional(ds, vars_dict, dimension_coord_names)
                if is_rotated_longitude_latitude(ds, vars_dict, attrs_dict, dimname, grid_mapping_key)
                    r = _roated_longitude_latudtude_dim(ds, vars_dict, attrs_dict, dimension_coord_names, dimname, crs)
                    return r
                elseif is_unalligned(ds, vars_dict, attrs_dict, dimension_coord_names, dimname, grid_mapping_key)
                    u = _unalligned_dim(ds, vars_dict, attrs_dict, dimension_coord_names, dimname, crs)
                    return u
                else
                    return _standard_dim(ds, vars_dict, attrs_dict, dimname, crs)
                end
            else
                return _merged_dim(ds, vars_dict, attrs_dict, dimension_coord_names, dimname, crs)
            end
        end
    end
    # TODO: allow Categorical, Sampled etc to have a nested dimension
    # remaining_coord_names = setdiff(dimension_coord_names, [linked_coord_name])
    # subdim = _cdm_dim(ds, vars_dict, attrs_dict, geometry_dict, geometry_key, var_name, remaining_coord_names, dimname, crs)
end


function is_rotated_longitude_latitude(ds, vars_dict, attrs_dict, dimname, grid_mapping_key)
    if !isnothing(grid_mapping_key)
        grid_mapping = _get_attr!(attrs_dict, vars_dict, ds, grid_mapping_key)
        grid_mapping_name = grid_mapping["grid_mapping_name"]
        if grid_mapping_name == "rotated_latitude_longitude" && get(dim_attrs, "standard_name", "") in ("grid_latitude", "grid_longitude")
            return true
        end
    end
    return false
end
function is_unalligned(ds, vars_dict, attrs_dict, dimension_coord_names, dimname, grid_mapping_key)
    units = map(dimension_coord_names) do coord_name
        coord_attrs = _get_attr!(attrs_dict, vars_dict, ds, coord_name)
        get(coord_attrs, "units", "")
    end
    return "degrees_north" in units && "degrees_east" in units
end
# Dimension categorisation
is_char_categorical(dim_var) = ndims(dim_var) > 1 && eltype(dim_var) <: Char
is_climatology(dim_attrs::AbstractDict) = haskey(dim_attrs, "climatology")
function is_multidimensional(ds, vars_dict, dimension_coord_names) 
    dim_vars = map(dimension_coord_names) do d
        _get_var!(vars_dict, ds, d)
    end
    return any(v -> length(CDM.dimnames(v)) > 1, dim_vars)
end
function is_geometry(ds, vars_dict, attrs_dict, geometry_key, dimname)
    isnothing(geometry_key) && return false
    geometry_attrs = _get_attr!(attrs_dict, vars_dict, ds, geometry_key)
    # Check that the geometry is for this dimension
    is_geometry_dimension = if haskey(geometry_attrs, "node_count")
        # Multi-part geometries 
        dimname in CDM.dimnames(_get_var!(vars_dict, ds, geometry_attrs["node_count"]))
    else                                             
        # Points                                     
        dimname in CDM.dimnames(_get_var!(vars_dict, ds, first(_cf_node_coordinate_names(geometry_attrs))))
    end
end

# Lookup/Dimension generation
function _standard_dim(ds, vars_dict, attrs_dict, dimname, crs)
    dim_attrs = _get_attr!(attrs_dict, vars_dict, ds, dimname)
    if is_climatology(dim_attrs)
        return _climatology_dim(ds, vars_dict, attrs_dict, dim_attrs, dimname)
    else
        D = _cdm_dimtype(dim_attrs, dimname)
        lookup = _cdm_lookup(ds, dim_attrs, dimname, D, crs)
        return D(lookup)
    end
end

function _discrete_axis_dim(ds, dimname)
    # It must be a "discrete axis" (CF 4.5) so NoLookup
    len = CDM.dim(ds, dimname)
    lookup = NoLookup(Base.OneTo(len))
    D = _cdm_dimtype(NoMetadata(), dimname)
    return D(lookup)
end

function _geometry_dim(ds, vars_dict, attrs_dict, geometry_dict, geometry_key, crs)
    geometry_attrs = _get_attr!(attrs_dict, vars_dict, ds, geometry_key)
    geoms = get!(geometry_dict, geometry_key) do
        _cf_geometry_decode(ds, geometry_attrs)
    end
    lookup = GeometryLookup(geoms, (X(), Y()); crs)
    D = Geometry
    return D(lookup)
end

function _climatology_dim(ds, vars_dict, attrs_dict, dim_attrs, name)
    periodfuncs = (year, month, day, hour, minute, second)

    # Get needed variables
    dim_var = ds[name]
    climatology_bounds = collect(CDM.cfvariable(ds, dim_attrs["climatology"]))

    # Detect the cycling pattern from the difference between periods in a step
    step_diff = map(periodfuncs) do f
        start = f(climatology_bounds[1, 1])
        stop = f(climatology_bounds[2, 1])
        stop - start
    end
    cycle, duration = if step_diff[1] != 0 && any(!=(0), step_diff[3:6]) # all change so Year and Day both cycle
        (Year(1), Day(1) => Month(1)), Year(step_diff[1])
    elseif step_diff[1] != 0 # year changes but not day/hour/minute/second so Year is the cycle
        Year(1), Year(step_diff[1])
    elseif any(!=(0), step_diff[3:6]) # day/hour/minute/second changes so Day is the cycle
        Day(1), Day(step_diff[3])
    else
        # There is no full cycle, but set it to 1 year for read/write round trip of the climatology bounds. 
        # It will have no effect on selectors.
        Year(1), Year(step_diff[1])
    end

    # Define lookup traits and values
    data = ndims(dim_var) == 0 ? [dim_var[]] : collect(dim_var)
    order = ForwardOrdered() # CF default ?
    sampling = Intervals(Center()) # CF default
    bounds = climatology_bounds[begin], climatology_bounds[end]
    # Reuse the vector but subtract duration from upper bounds
    # This converts climatology bounds to normal bounds
    climatology_bounds[2, :] .-= duration
    span = Explicit(climatology_bounds)
    lookup = Cyclic(data; cycle, order, sampling, span, bounds)

    # Detect dimention
    D = _cdm_dimtype(dim_attrs, name)
    return D(lookup)
end

function _char_categorical_dim(ds, dim_var, dim_attrs, name)
    chars = collect(dim_var)
    # Join char dimension as String, stripping null terminators
    strings = String.(strip.(join.(eachslice(chars; dims=2)), '\0'))
    metadata = _metadatadict(sourcetrait(ds), dim_attrs)

    lookup = Categorical(strings; order=Unordered(), metadata)
    D = _cdm_dimtype(dim_attrs, name)
    return D(lookup)
end

# Must have "standard_name" in dimension attributes
function _rotated_longitude_latudtude_dim(ds, vars_dict, attrs_dict, dimension_coord_names, dimname, crs)
    lookup = _unaligned_lookup(ds, vars_dict, attrs_dict, dimension_coord_names, dimname, crs)
    standard_name = _get_attr!(attrs_dict, vars_dict, ds)["standard_name"]
    D = if standard_name == "grid_latitude"
        Y
    elseif standard_name == "grid_longitude"
        X
    else
        error("standard_name $standard_name not recognized")
    end
    return D(lookup)
end
# Must have "degrees_north" and "degrees_east" in coordinate units
function _unalligned_dim(ds, vars_dict, attrs_dict, dimension_coord_names, dimname, crs)
    dim_attrs = _get_attr!(attrs_dict, vars_dict, ds, dimname)
    D = _cdm_dimtype(dim_attrs, dimname)
    lookup = _unaligned_lookup(ds, vars_dict, attrs_dict, dimension_coord_names, dimname, D, crs)
    return D(lookup)
end

function _unaligned_lookup(ds, vars_dict, attrs_dict, dimension_coord_names, dimname, D, crs)
    # Rotations. For now we don't use ArrayLookup
    lon_key = lat_key = ""
    for coord_name in dimension_coord_names
        coord_attrs = _get_attr!(attrs_dict, vars_dict, ds, coord_name)
        units = get(coord_attrs, "units", nothing)
        if units == "degrees_north"
            lat_key = coord_name
        elseif units == "degrees_east"
            lon_key = coord_name
        end
    end
    lat_metadata, lon_metadata = map((lat_key, lon_key)) do key
        _metadatadict(sourcetrait(ds), _get_attr!(attrs_dict, vars_dict, ds, key))
    end
    dims = Lat(AutoValues(); metadata=lat_metadata), 
           Lon(AutoValues(); metadata=lon_metadata)
    data = if haskey(ds, dimname)
        _get_attr!(attrs_dict, vars_dict, ds, dimname)
    else
        Base.OneTo(CDM.dim(ds, dimname))
    end
    if D <: X
        matrix = collect(_get_var!(vars_dict, ds, lon_key))
        return ArrayLookup(matrix; data, dims)
    elseif D == Y
        matrix = collect(_get_var!(vars_dict, ds, lat_key))
        return ArrayLookup(matrix; data, dims)
    end
end

function _merged_dim(ds, vars_dict, attrs_dict, dimension_coord_names, dimname, crs)
    # Generate a MergedLookup for multi-coordinate dimensions
    dims = map(dimension_coord_names) do d
        dim_attrs = _get_attr!(attrs_dict, vars_dict, ds, d)
        D = _cdm_dimtype(dim_attrs, d)
        lookup = _cdm_lookup(ds, dim_attrs, d, D, crs)
        D(lookup)
    end
    D = _cdm_dimtype(NoMetadata(), dimname)
    merged_data = collect(zip(DD.lookup.(dims)...))
    # Use a MergedLookup with all dimensions
    return D(MergedLookup(merged_data, Tuple(dims)))
end

# TODO don't load all keys here with _layers
_name_or_firstname(ds::AbstractDataset, name) = Symbol(name)
function _name_or_firstname(ds::AbstractDataset, name::Union{Nothing,NoKW})
    names = _layer_names(ds)
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
function _cdm_lookup(ds::AbstractDataset, attr, dimname, D::Type, crs)
    var = ds[dimname]
    data = Missings.disallowmissing(collect(var))
    # Length one variables may be zero dimensional
    data = ndims(data) == 0 ? reshape(data, 1) : data
    metadata = _metadatadict(sourcetrait(ds), attr)
    _cdm_lookup(data, ds, var, attr, dimname, D, metadata, crs)
end
# For unknown types we just make a Categorical lookup
function _cdm_lookup(data::AbstractArray, ds::AbstractDataset, var, attr, dimname, D::Type, metadata, crs)
    Categorical(data; order=Unordered(), metadata=metadata)
end
# For Number and AbstractTime we generate order/span/sampling
# We need to include `Missing` in unions in case `_FillValue` is used
# on coordinate variables in a file and propagates here.
function _cdm_lookup(
    data::AbstractArray{<:Union{Missing,Number,Dates.AbstractTime}},
    ds::AbstractDataset, var, attr, dimname, D::Type, metadata, crs
)
    # Assume the locus is at the center of the cell if boundaries aren't provided.
    # http://cfconventions.org/cf-conventions/cf-conventions.html#cell-boundaries
    order = LA.orderof(data)
    # Detect lat/lon
    span, sampling = if eltype(data) <: Union{Missing,Dates.AbstractTime}
        _cdm_period(data, metadata)
    else
        _cdm_span(data, order)
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
            locus = if mapreduce(==, &, view(boundsmatrix, 1, :), data)
                Start()
            elseif mapreduce(==, &, view(boundsmatrix, 2, :), data)
                End()
            else
                # TODO better checks here
                Center()
            end
            Explicit(boundsmatrix), Intervals(locus)
        end
    end
    return _cdm_lookup(data, D, order, span, sampling, metadata, crs)
end
# For X and Y use a Mapped <: AbstractSampled lookup
function _cdm_lookup(
    data, D::Type{<:Union{<:XDim,<:YDim}}, order::Order, span, sampling, metadata, crs; dims=nothing
)
    # If units are degrees north/east, use EPSG:4326 as the mapped crs
    units = get(metadata, "units", "")
    mappedcrs = if (units in CF_DEGREES_NORTH || units in CF_DEGREES_EAST)
        EPSG(4326)
    end
    # Additionally, crs and mappedcrs should be identical for Regular lookups
    if isnokwornothing(crs) && span isa Regular
        crs = mappedcrs
    end
    dim = DD.basetypeof(D)()
    return Mapped(data; order, span, sampling, metadata, crs, mappedcrs, dim, dims)
end
# Band dims have a Categorical lookup, with order
# This is not CF, just for consistency with GDAL
function _cdm_lookup(data, D::Type{<:Band}, order::Order, span, sampling, metadata, crs)
    Categorical(data, order, metadata)
end
# Otherwise use a Sampled lookup
function _cdm_lookup(data, D::Type, order::Order, span, sampling, metadata, crs)
    Sampled(data, order, span, sampling, metadata)
end

function _cdm_span(data, order)
    # Handle a length 1 data
    length(data) == 1 && return Regular(zero(eltype(data))), Points()

    step = if eltype(data) <: AbstractFloat
        # Calculate step, avoiding as many floating point errors as possible
        st = Base.step(Base.range(Float64(first(data)), Float64(last(data)); length = length(data)))
        st_rd = round(st, digits = Base.floor(Int,-log10(eps(eltype(data))))) # round to nearest digit within machine epsilon
        isapprox(st_rd, st; atol = eps(eltype(data))) ? st_rd : st # keep the rounded number if it is very close to the original
    else
        data[2] - data[1]
    end
    for i in 2:length(data)-1
        # If any step sizes don't match, its Irregular
        if !(data[i+1] - data[i] ≈ step)
            bounds = if length(data) > 1
                beginhalfcell = abs((data[2] - data[1]) * 0.5)
                endhalfcell = abs((data[end] - data[end-1]) * 0.5)
                if LA.isrev(order)
                    data[end] - endhalfcell, data[1] + beginhalfcell
                else
                    data[1] - beginhalfcell, data[end] + endhalfcell
                end
            else
                data[1], data[1]
            end
            return Irregular(bounds), Points()
        end
    end
    # Otherwise regular
    return Regular(step), Points()
end

# delta_t and ave_period are not CF standards, but CDC
function _cdm_period(data, metadata::Metadata{<:CDMsource})
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

_cdm_shiftlocus(dim::Dimension) = _cdm_shiftlocus(lookup(dim), dim)
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

    varname = if isnokw(name) || string(name) == ""
        UNNAMED_FILE_KEY
    else
        string(name)
    end

    if !isempty(dims(A, x -> lookup(x) isa GeometryLookup))
        attrib["geometry"] = "geometry_container"
    end
    dimnames = map(_cf_name, dims(A))
    var = CDM.defVar(ds, varname, eltype, dimnames; attrib, chunksizes, kw...)

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
    dimname = _cf_name(dim)
    haskey(ds.dim, dimname) && return nothing
    CDM.defDim(ds, dimname, length(dim))
    _def_lookup_var!(ds, dim, dimname)
    return nothing
end 
_def_lookup_var!(ds::AbstractDataset, dim::Dimension{<:AbstractNoLookup}, dimname) = nothing
function _def_lookup_var!(ds::AbstractDataset, dim::Dimension{<:ArrayLookup}, dimname)  
    # TODO what do we call the geometry variable, what if more than 1
    return nothing
end
function _def_lookup_var!(ds::AbstractDataset, dim::Dimension{<:GeometryLookup}, dimname)  
    l = lookup(dim)
    geom = _cf_geometry_encode(lookup(dim))
    # Define base geometry attributes
    coord_names_raw = map(dims(l)) do d
        _cf_name(d)
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
        lon_varname = "lon"
        lat_varname = "lat"
        node_count_varname = dimname * "_node_count"
        CDM.defDim(ds, node_dimname, length(geom.x))
        CDM.defVar(ds, node_count_varname, geom.node_count, (dimname,))
        CDM.defVar(ds, lon_varname, geom.lon, (dimname,); 
            attrib=["standard_name" => "longitude", "nodes" => coord_names[1]]
        )
        CDM.defVar(ds, lat_varname, geom.lat, (dimname,); 
            attrib=["standard_name" => "latitude", "nodes" => coord_names[2]]
        )
        attrib["node_count"] = node_count_varname
        attrib["coordinates"] = "lat lon"
    else
        # Node count is dimension length
        node_dimname = dimname
    end
    if haskey(geom, :part_node_count)
        # We also have geometry parts
        part_dimname = dimname * "_part"
        part_varname = dimname * "_part_node_count"
        CDM.defDim(ds, part_dimname, length(geom.part_node_count))
        CDM.defVar(ds, part_varname, geom.part_node_count, (part_dimname,))
        attrib["part_node_count"] = part_varname
    end
    if haskey(geom, :interior_ring)
        # And interior rings
        ring_varname = dimname * "_interior_ring"
        CDM.defVar(ds, ring_varname, geom.interior_ring, (part_dimname,))
        attrib["interior_ring"] = ring_varname
    end
    # Define coordinate variables
    CDM.defVar(ds, coord_names[1], geom.x, (node_dimname,))
    CDM.defVar(ds, coord_names[2], geom.y, (node_dimname,))
    # TODO: add z and m, if present
    # And the geometry container
    CDM.defVar(ds, "geometry_container", fill(0), (); attrib)
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
    CDM.defVar(ds, dimname, collect(lookup(dim)), (dimname,); attrib)
    return nothing
end
    
_cf_crs(ds, grid_mapping_key::Nothing) = nothing
function _cf_crs(ds, grid_mapping_key::String)
    # TODO handle multi-crs
    if occursin(' ', grid_mapping_key)
        _split_attribute(grid_mapping_key)[1][1]
    end
    haskey(ds, grid_mapping_key) || return nothing
    grid_mapping_var = CDM.variable(ds, grid_mapping_key)
    grid_mapping_attrib = CDM.attribs(grid_mapping_var)
    if haskey(grid_mapping_attrib, "crs_wkt")
        WellKnownText2(GeoFormatTypes.CRS(), grid_mapping_attrib["crs_wkt"])
    elseif haskey(grid_mapping_attrib, "spatial_epsg")
        EPSG(parse(Int, grid_mapping_attrib["spatial_epsg"]))
    elseif haskey(grid_mapping_attrib, "proj4string")
        ProjString(grid_mapping_attrib["proj4string"])
    else
        if haskey(grid_mapping_attrib, "grid_mapping_name")
            grid_mapping_name = grid_mapping_attrib["grid_mapping_name"]
            if grid_mapping_name == "latitude_longitude"
                return EPSG(4326)
            end
        end
        nothing # unparseable - TODO get CFProjections.jl to do it.
    end
end

function _split_attribute(attribute::String)
    items = split(attribute, ' ')
    starts = findall(contains(':'), items)
    stops = append!(starts[2:end] .- 1, length(items))
    map(starts, stops) do s, p
        items[s][1:end-1] => items[s+1:p]
    end
end

function _grid_mapping_keys(grid_mapping_key::String)
    if occursin(' ', grid_mapping_key)
        first.(_split_attribute(grid_mapping_key))
    else
        [grid_mapping_key]
    end
end

# Get a variable and store it in a dict
_get_var!(dict::AbstractDict, ds, key::String) = get!(() -> CDM.variable(ds, key), dict, key)

# Get atributes and store them in a dict
function _get_attr!(attr_dict::AbstractDict, var_dict::AbstractDict, ds, key::String)
    get!(attr_dict, key) do 
        var = _get_var!(var_dict, ds, key)
        CDM.attribs(var)
    end
end

# We need a custom `eltype` method because Char is convered to String
_eltype(::CDMsource, var::AbstractArray{Char}) = String
_eltype(::Source, var) = eltype(var)

# TODO keep just the core of this function here
function _raster(ds::CDM.AbstractDataset;
    dims=nokw,
    refdims=nokw,
    name=nokw,
    group=nokw,
    filename=filename(ds),
    metadata=nokw,
    missingval=nokw,
    crs=nokw,
    mappedcrs=nokw,
    source=sourcetrait(ds),
    replace_missing=nokw,
    coerce=convert,
    scaled::Union{Bool,NoKW}=nokw,
    verbose::Bool=true,
    write::Bool=false,
    lazy::Bool=false,
    dropband::Bool=true,
    checkmem::Bool=CHECKMEM[],
    raw::Bool=false,
    mod=nokw,
    f=identity,
)::Raster
    _maybe_warn_replace_missing(replace_missing)
    # `raw` option will ignore `scaled` and `missingval`
    scaled, missingval = _raw_check(raw, scaled, missingval, verbose)
    # TODO use a clearer name than filekey
    name1 = filekey(ds, name)
    (; names_vec, layers_vec, layerdims_vec, layermetadata_vec, dim_dict, refdim_dict) = 
        _organise_dataset(ds, [string(name1)], group)
    var = layers_vec[1]
    # Detect the source from filename
    # Open the dataset and variable specified by `name`, at `group` level if provided
    # At this level we do not apply `mod`.
    refdims = isnokw(refdims) ? Tuple(values(refdim_dict)) : refdims
    metadata_out = isnokw(metadata) ? Dict(layermetadata_vec[1]) : metadata
    missingval_out = _read_missingval_pair(var, metadata_out, missingval)
    pop!(metadata_out, "_FillValue", nothing)
    pop!(metadata_out, "missing_value", nothing)
    # Generate mod for scaling
    mod = isnokw(mod) ? _mod(_eltype(source, var), metadata_out, missingval_out; scaled, coerce) : mod
    # Define or load the data array
    data_out = if lazy
        # Define a lay FileArray
        FileArray{typeof(source)}(var, filename;
            name=name1, group, mod, write
        )
    else
        modvar = _maybe_modify(var, mod)
        # Check the data will fit in memory
        checkmem && _checkobjmem(modvar)
        # Move the modified array to memory
        Array(modvar)
    end
    # Generate dims
    dims_out = isnokw(dims) ? DD.dims(Tuple(values(dim_dict)), layerdims_vec[1]) : dims
    # Return the data to the parent function
    mv_outer = _outer_missingval(mod)
    # Use name or an empty Symbol
    name_out = name1 isa Union{NoKW,Nothing} ? Symbol("") : Symbol(name1)
    # Define the raster
    raster = Raster(data_out, dims_out, refdims, name_out, metadata_out, missingval_out)
    # Maybe drop a single band dimension
    return _maybe_drop_single_band(raster, dropband, lazy)
end
