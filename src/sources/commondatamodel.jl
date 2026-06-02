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
DiskArrays.eachchunk(A::DiskCharToString) = DiskArrays.GridChunks(eachchunk(parent(A)).chunks[2:end])
function DiskArrays.readblock!(A::DiskCharToString, dest, I...)
    src = parent(A)[:, I...]
    for I in CartesianIndices(dest)
        dest[I] = _cf_chars_to_string(view(src, :, I))
    end
    return dest
end

_cf_chars_to_string(chars) = String(strip(join(chars), '\0'))

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
    quantization_names = String[]
    climatology_names = String[]
    # Sometimes coordinates are not included in var attributes
    # So we check standard names to know vars are lookups and not data
    standard_names = String[]
    for n in var_names
        var = _get_var!(vars_dict, ds, n)
        attr = _get_attr!(attrs_dict, vars_dict, ds, n)
        # Coordinate variables. UGRID `node_coordinates` (CF 5.21 mesh
        # topology) are not filtered here on purpose: they are 1D coord vars
        # the writer needs to put back. They surface as layers, which the
        # mesh's `node_coordinates` attribute then references correctly.
        haskey(attr, "coordinates") && append!(coordinate_names, collect(_cdm_coordinates(attr)))
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
        standard_names
    )
end

_cdm_coordinates(attr) = 
    haskey(attr, "coordinates") ? _split_attribute(attr["coordinates"]) : String[]
_cdm_node_coordinates(attr) =
    haskey(attr, "node_coordinates") ? _split_attribute(attr["node_coordinates"]) : String[]

_split_attribute(attrib) = 
    isempty(attrib) ? String[] : String.(split(attrib, ' '))

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
    dim_dict = OrderedDict{String,Dimension}()
    # Refdims need a consistent order as it is not changed later
    refdim_dict = OrderedDict{String,Dimension}()
    # Get the layer names to target
    layer_names = isnokw(names) ? _layer_names(ds, var_names, vars_dict, attrs_dict) : names
    layer_dimnames_vec = Vector{Tuple}(undef, length(layer_names))

    # Define output data and metadata vars
    output_layers_vec = Vector{AbstractArray}(undef, length(layer_names))
    output_layerdims_vec = Vector{Tuple}(undef, length(layer_names))
    output_attrs_vec = Vector{Any}(undef, length(layer_names))
    used_layers = trues(length(layer_names))

    # Loop over the layers we want to load as rasters
    # As we go, we add dimensions to dim_dict to be used in other layers
    for (i, var_name) in enumerate(layer_names)
        # Get the variable and its attributes
        var = get(() -> CDM.variable(ds, var_name), vars_dict, var_name)
        attr = get(() -> CDM.attribs(var), attrs_dict, var_name)
        isdomain = false
        if eltype(var) <: Char && ndims(var) == 0
            if haskey(attr, "dimensions")
                layer_dimnames_vec[i] = dimnames = Tuple(_split_attribute(attr["dimensions"]))
                # Remove this entry from layer vecs
                isdomain = true
            end
        else
            # Get its dimensions
            layer_dimnames_vec[i] = dimnames = CDM.dimnames(var)
        end
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
        layerdims = map(dimnames) do dimname
            get!(dim_dict, dimname) do
                # Resolve crs per-dim. For multi-grid_mapping files (CF 5.10
                # has "crsOSGB: x y crsWGS84: lat lon"), the loader earlier
                # returned a Vector{Pair{CRS, Vector{String}}} - which made
                # crs(stack) a vector, breaking everything downstream that
                # expects a single GeoFormat. Per-dim resolution gives each
                # dim's lookup the single CRS that applies to its coord var.
                dim_crs, dim_mappedcrs = _resolve_dim_crs(crs, dimname, coord_names)
                _cdm_dim(
                    ds, vars_dict, attrs_dict, geometry_dict, geometry_key, var_name, coord_names, dimname, dim_crs, dim_mappedcrs, grid_mapping_key
                )
            end
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

        # Format dimensions
        unformatted_dims = map(dimnames) do dimname
            dim_dict[dimname]
        end
        if isdomain
            formatted_dims = format(unformatted_dims)
            used_layers[i] = false
        else
            output_layerdims_vec[i] = layerdims
            maybewrappedvar = if eltype(var) <: Char && length(layerdims) == (ndims(var) - 1)
                DiskCharToString(var)
            else
                var
            end
            formatted_dims = format(unformatted_dims, maybewrappedvar)
            # Finalise output variables and attributes
            output_layers_vec[i] = maybewrappedvar
            output_attrs_vec[i] = attr
        end
        map(dimnames, formatted_dims) do dimname, d
            dim_dict[dimname] = d
        end
    end

    return (; 
        names_vec=layer_names[used_layers], 
        layers_vec=output_layers_vec[used_layers], 
        layerdims_vec=output_layerdims_vec[used_layers], 
        layermetadata_vec=output_attrs_vec[used_layers], 
        dim_dict,
        refdim_dict,
    )
end
_organise_dataset(ds::AbstractDataset, names, group) =
    _organise_dataset(ds.group[group], names, nokw)

# Simple version for reading dimensions from a single variable (not full dataset analysis)
function _cdm_dim(ds, dimname, crs, mappedcrs)
    if haskey(ds, dimname)
        # There's a coordinate variable matching the dimension name
        vars_dict = Dict{String,Any}()
        attrs_dict = Dict{String,Any}()
        return _standard_dim(ds, vars_dict, attrs_dict, dimname, crs)
    else
        # No coordinate variable - use discrete axis (NoLookup)
        return _discrete_axis_dim(ds, dimname)
    end
end

function _cdm_dim(ds, vars_dict, attrs_dict, geometry_dict, geometry_key, var_name, coord_names, dimname, crs, mappedcrs, grid_mapping_key)
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
            if haskey(ds, dimname) && dimname ∉ dimension_coord_names
                # There's a dimension variable AND an auxiliary coordinate - create MergedLookup
                # CF 6.2: dimension variable (sigma) + auxiliary coordinate (model_level)
                all_coord_names = [dimname; dimension_coord_names]
                return _merged_dim(ds, vars_dict, attrs_dict, all_coord_names, dimname, crs)
            elseif haskey(ds, dimname)
                # Only one coordinate, so just use it as the dimension/lookup
                return _standard_dim(ds, vars_dict, attrs_dict, dimname, crs)
            else
                return _standard_dim(ds, vars_dict, attrs_dict, dimension_coord_names[1], crs)
            end
        else
            first_coord_name = dimension_coord_names[1]
            first_coord_var = _get_var!(vars_dict, ds, first_coord_name)
            # Handle multiple char categorical dims (CF 6.1.2: taxon_lsid + taxon_name)
            if is_char_categorical(first_coord_var)
                return _char_categorical_merged_dim(ds, vars_dict, attrs_dict, dimension_coord_names, dimname)
            end
            if is_multidimensional(ds, vars_dict, dimension_coord_names)
                if is_rotated_longitude_latitude(ds, vars_dict, attrs_dict, dimname, grid_mapping_key)
                    r = _rotated_longitude_latudtude_dim(ds, vars_dict, attrs_dict, dimension_coord_names, dimname, crs, mappedcrs)
                    return r
                elseif is_unalligned(ds, vars_dict, attrs_dict, dimension_coord_names, dimname, grid_mapping_key)
                    u = _unalligned_dim(ds, vars_dict, attrs_dict, dimension_coord_names, dimname, crs, mappedcrs)
                    return u
                else
                    return _standard_dim(ds, vars_dict, attrs_dict, dimname, crs)
                end
            else
                if haskey(ds, dimname)
                    # This is something we're not handling properly yet, try to warn
                    dim_attrs = _get_attr!(attrs_dict, vars_dict, ds, dimname)
                    _check_formula_terms(dim_attrs)
                end
                # Only one coordinate, so just use it as the dimension/lookup
                return _merged_dim(ds, vars_dict, attrs_dict, dimension_coord_names, dimname, crs)
            end
        end
    end
    # TODO: allow Categorical, Sampled etc to have a nested dimension
    # remaining_coord_names = setdiff(dimension_coord_names, [linked_coord_name])
    # subdim = _cdm_dim(ds, vars_dict, attrs_dict, geometry_dict, geometry_key, var_name, remaining_coord_names, dimname, crs)
end

# Multi-grid_mapping resolution. CF 5.10 has
#   temp:grid_mapping = "crsOSGB: x y crsWGS84: lat lon"
# which the loader turns into a Vector{Pair{CRS, Vector{String}}}. Each dim
# only wants the single CRS that applies to its coord var; the OTHER CRS in
# the multi-mapping (here WGS84, scoped to aux lat/lon) becomes the
# `mappedcrs` of an unaligned/ProjectedArrayLookup dim.
_resolve_dim_crs(crs, dimname, aux_coord_names) = (crs, nothing)
function _resolve_dim_crs(crs::AbstractVector, dimname, aux_coord_names)
    primary = nothing
    mapped = nothing
    for (one_crs, coord_names) in crs
        if dimname in coord_names
            primary = one_crs
        elseif any(c -> c in coord_names, aux_coord_names)
            mapped = one_crs
        end
    end
    return primary, mapped
end

function _check_formula_terms(dim_attrs)
    haskey(dim_attrs, "formula_terms") && @warn """
    formula_terms are not yet recognised or used by Rasters.jl. 
    This may change at any time in future without a breaking version."
    """
end

function is_rotated_longitude_latitude(ds, vars_dict, attrs_dict, dimname, grid_mapping_key)
    isnothing(grid_mapping_key) && return false
    # grid_mapping_key may be a complex string with multiple mappings, check if it's a valid variable
    haskey(ds, grid_mapping_key) || return false
    grid_mapping = _get_attr!(attrs_dict, vars_dict, ds, grid_mapping_key)
    grid_mapping_name = get(grid_mapping, "grid_mapping_name", nothing)
    isnothing(grid_mapping_name) && return false
    if grid_mapping_name == "rotated_latitude_longitude"
        dim_attrs = _get_attr!(attrs_dict, vars_dict, ds, dimname)
        if get(dim_attrs, "standard_name", "") in ("grid_latitude", "grid_longitude")
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

# Check if 2D coordinate variables have polygon bounds (CF 7.2)
function has_polygon_bounds(ds, vars_dict, attrs_dict, dimension_coord_names)
    for coord_name in dimension_coord_names
        coord_attrs = _get_attr!(attrs_dict, vars_dict, ds, coord_name)
        if haskey(coord_attrs, "bounds")
            bounds_name = coord_attrs["bounds"]
            if haskey(ds, bounds_name)
                bounds_var = _get_var!(vars_dict, ds, bounds_name)
                bounds_dimnames = CDM.dimnames(bounds_var)
                # Bounds should have 3 dims: the coordinate dims plus a vertex dim
                # with at least 3 vertices for polygons
                if length(bounds_dimnames) >= 3
                    # Find the vertex dimension (not in coordinate dimnames)
                    coord_var = _get_var!(vars_dict, ds, coord_name)
                    coord_dimnames = CDM.dimnames(coord_var)
                    vertex_dim = first(filter(d -> d ∉ coord_dimnames, bounds_dimnames))
                    nv = CDM.dim(ds, vertex_dim)
                    if nv >= 3
                        return true
                    end
                end
            end
        end
    end
    return false
end

# Create polygon geometries from 2D bounds arrays (CF 7.2)
function _create_polygons_from_bounds(ds, vars_dict, attrs_dict, lat_key, lon_key; crs=nothing)
    lat_attrs = _get_attr!(attrs_dict, vars_dict, ds, lat_key)
    lon_attrs = _get_attr!(attrs_dict, vars_dict, ds, lon_key)

    lat_bounds_name = lat_attrs["bounds"]
    lon_bounds_name = lon_attrs["bounds"]

    lat_bnds = collect(_get_var!(vars_dict, ds, lat_bounds_name))
    lon_bnds = collect(_get_var!(vars_dict, ds, lon_bounds_name))

    # lat_bnds and lon_bnds are 3D: (nv, dim1, dim2)
    nv, ni, nj = size(lat_bnds)

    # Create a polygon for each cell
    # Use the same approach as _cf_geometry_decode for polygons
    polygons = Vector{Any}(undef, ni * nj)
    idx = 0
    for j in 1:nj, i in 1:ni
        idx += 1
        # Get the vertices for this cell in (lon, lat) order
        vertices = [(lon_bnds[v, i, j], lat_bnds[v, i, j]) for v in 1:nv]
        # Close the polygon by adding the first vertex at the end
        push!(vertices, vertices[1])
        # Create the polygon
        ring = GI.LinearRing(vertices; crs)
        polygons[idx] = GI.Polygon([ring]; crs)
    end

    return polygons, (ni, nj)
end

# Lookup/Dimension generation
function _standard_dim(ds, vars_dict, attrs_dict, dimname, crs)
    dim_attrs = _get_attr!(attrs_dict, vars_dict, ds, dimname)
    _check_formula_terms(dim_attrs)
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
    # Climatology bounds inherit time units from the parent dimension variable
    # per CF convention. cfvariable won't auto-decode without a `units` attribute
    # on the bounds variable itself, so pass the parent's units explicitly.
    bounds_kw = if haskey(dim_attrs, "units")
        (; units=dim_attrs["units"])
    else
        (;)
    end
    climatology_bounds = collect(CDM.cfvariable(ds, dim_attrs["climatology"]; bounds_kw...))

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

# Handle multiple char categorical coordinates (CF 6.1.2: taxon_lsid + taxon_name)
function _char_categorical_merged_dim(ds, vars_dict, attrs_dict, dimension_coord_names, dimname)
    # Create a lookup for each char categorical coordinate
    dims = map(dimension_coord_names) do coord_name
        coord_var = _get_var!(vars_dict, ds, coord_name)
        coord_attrs = _get_attr!(attrs_dict, vars_dict, ds, coord_name)
        chars = collect(coord_var)
        # Join char dimension as String, stripping null terminators
        strings = String.(strip.(join.(eachslice(chars; dims=2)), '\0'))
        metadata = _metadatadict(sourcetrait(ds), coord_attrs)
        lookup = Categorical(strings; order=Unordered(), metadata)
        D = _cdm_dimtype(coord_attrs, coord_name)
        D(lookup)
    end
    # Create MergedLookup with all dimensions
    D = _cdm_dimtype(NoMetadata(), dimname)
    merged_data = collect(zip(DD.lookup.(dims)...))
    return D(MergedLookup(merged_data, Tuple(dims)))
end

# Must have "standard_name" in dimension attributes
function _rotated_longitude_latudtude_dim(ds, vars_dict, attrs_dict, dimension_coord_names, dimname, crs, mappedcrs=nothing)
    standard_name = _get_attr!(attrs_dict, vars_dict, ds, dimname)["standard_name"]
    D = if standard_name == "grid_latitude"
        Y
    elseif standard_name == "grid_longitude"
        X
    else
        error("standard_name $standard_name not recognized")
    end
    lookup = _unaligned_lookup(ds, vars_dict, attrs_dict, dimension_coord_names, dimname, D, crs, mappedcrs)
    return D(lookup)
end
# Must have "degrees_north" and "degrees_east" in coordinate units
function _unalligned_dim(ds, vars_dict, attrs_dict, dimension_coord_names, dimname, crs, mappedcrs=nothing)
    # dimname may be just a dimension without a corresponding variable
    dim_attrs = haskey(ds, dimname) ? _get_attr!(attrs_dict, vars_dict, ds, dimname) : NoMetadata()
    D = _cdm_dimtype(dim_attrs, dimname)
    lookup = _unaligned_lookup(ds, vars_dict, attrs_dict, dimension_coord_names, dimname, D, crs, mappedcrs)
    return D(lookup)
end

function _unaligned_lookup(ds, vars_dict, attrs_dict, dimension_coord_names, dimname, D, crs, mappedcrs=nothing)
    # Identify the lat / lon aux coord vars by CF metadata, not by name.
    # `units = "degrees_east"` (or `standard_name = "longitude"`) -> lon;
    # `units = "degrees_north"` (or `standard_name = "latitude"`)  -> lat.
    lon_key = lat_key = ""
    for coord_name in dimension_coord_names
        coord_attrs = _get_attr!(attrs_dict, vars_dict, ds, coord_name)
        units = get(coord_attrs, "units", nothing)
        std = get(coord_attrs, "standard_name", nothing)
        if units == "degrees_north" || std == "latitude"
            lat_key = coord_name
        elseif units == "degrees_east" || std == "longitude"
            lon_key = coord_name
        end
    end
    lat_metadata, lon_metadata = map((lat_key, lon_key)) do key
        _metadatadict(sourcetrait(ds), _get_attr!(attrs_dict, vars_dict, ds, key))
    end
    dims = Lat(AutoValues(); metadata=lat_metadata),
           Lon(AutoValues(); metadata=lon_metadata)
    if haskey(ds, dimname)
        dim = _cdm_dimtype(ds, dimname)()
        data = collect(_get_var!(vars_dict, ds, dimname))
    else
        dim = AutoDim()
        data = Base.OneTo(CDM.dim(ds, dimname))
    end
    # Decide which aux coord this lookup's matrix represents.
    # CF convention: an X-type outer dim pairs with the longitude aux coord
    # (degrees_east); a Y-type outer dim pairs with the latitude aux coord.
    # For generic `Dim`-typed outer dims (CF 7.2) there is no CF rule, so we
    # default to lon for the first such dim encountered and lat for the
    # second - the loader sees both aux coord vars and the writer reads
    # `aux_name` to pick the right var name back. This still loses one of
    # the two matrices per file (only the chosen aux coord is stored in
    # `matrix`); a richer ProjectedArrayLookup that stores both matrices
    # would address that.
    aux_name, matrix_key, paired_key = if D <: XDim
        :lon, lon_key, lat_key
    elseif D <: YDim
        :lat, lat_key, lon_key
    else
        # Plain Dim outer axis - no CF rule to pick. Default to lat as the
        # primary so the X/Y default (X=lon, Y=lat) stays consistent.
        :lat, lat_key, lon_key
    end
    matrix = collect(_get_var!(vars_dict, ds, matrix_key))
    # Load the paired aux coord matrix too (when present) so the writer can
    # put both back on disk - the on-disk file has both, but the in-memory
    # lookup only owned one before.
    paired_matrix = isempty(paired_key) ? nothing :
        collect(_get_var!(vars_dict, ds, paired_key))

    # Check for polygon bounds (CF 7.2) and create GeometryLookup if available
    geom_lookup = nothing
    if has_polygon_bounds(ds, vars_dict, attrs_dict, dimension_coord_names)
        polygons, _ = _create_polygons_from_bounds(ds, vars_dict, attrs_dict, lat_key, lon_key; crs)
        geom_lookup = GeometryLookup(polygons, (X(), Y()); crs)
    end

    return ProjectedArrayLookup(matrix; data, dims, crs, mappedcrs, geom_lookup, aux_name, paired_matrix)
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
function _cdm_lookup(data::AbstractVector{<:Union{AbstractString,Char}}, ds::AbstractDataset, var, attr, dimname, D::Type, metadata, crs)
    Categorical(data; order=Unordered(), metadata=metadata)
end
function _cdm_lookup(data::AbstractMatrix{Char}, ds::AbstractDataset, var, attr, dimname, D::Type, metadata, crs)
    Categorical(_cf_chars_to_string.(eachslice(data; dims=2)); order=Unordered(), metadata=metadata)
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
    # Replace raw float vectors with a `StableRange` when the span is regular,
    # so slicing/mosaicking/`Contains` stay bit-stable downstream. Done after
    # the bounds-locus check above, which relies on disk values matching bounds
    # bit-exactly via `==`.
    if span isa Regular && eltype(data) <: AbstractFloat && length(data) > 1
        data = StableRange(; start=first(data), step=val(span), length=length(data))
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
    if isnothing(mtch)
        return nothing
    else
        if length(mtch.captures) == 6
            vals = ntuple(Val{6}()) do i
                x = mtch.captures[i]
                # TODO can it actually be nothing?
                isnothing(x) ? 0 : parse(Int, x)
            end
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
# `get!` rather than `[]=` so we never overwrite a `standard_name` that came
# from the source file (e.g. CF 5.10's Z dim has standard_name=
# "height_above_reference_ellipsoid", which should survive a round-trip).
_cdm_set_axis_attrib!(atr, dim::X) = get!(atr, "axis", "X")
_cdm_set_axis_attrib!(atr, dim::Y) = get!(atr, "axis", "Y")
_cdm_set_axis_attrib!(atr, dim::Z) = (get!(atr, "axis", "Z"); get!(atr, "standard_name", "depth"))
_cdm_set_axis_attrib!(atr, dim::Ti) = (get!(atr, "axis", "T"); get!(atr, "standard_name", "time"))
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
        # Write the stack's full set of dimensions up-front. Without this,
        # `writevar!` is the only thing that calls `_def_dim_var!`, so any
        # dim that no layer uses gets dropped on round-trip - a "domain"
        # stack (CF 5.15/5.16/5.19, layers=()) loses all dims, and 5.17
        # loses its `Ti` dim because only `cell_area(cell)` references the
        # layer-bound dims.
        stack_aux = _def_dim_var!(ds, s)
        mods = map(keys(s)) do k
            writevar!(ds, source, s[k]; missingval=missingval[k], kw...)
        end
        # If any of the stack's dims are not covered by layer variables, write
        # a CF "domain" marker variable (CF section 5.7) so the loader knows
        # those dims belong to the stack. A scalar char with a `dimensions`
        # attribute is the CF convention.
        _maybe_write_domain!(ds, source, s, stack_aux)
        f(OpenStack{Source,K,T}(ds, mods))
    finally
        close(ds)
    end
    return filename
end

function _maybe_write_domain!(ds::AbstractDataset, source::CDMsource, s::AbstractRasterStack, stack_aux)
    stack_dims = dims(s)
    isempty(stack_dims) && return nothing
    layer_keys = keys(s)
    covered_names = if isempty(layer_keys)
        Set{Symbol}()
    else
        Set(DD.name(d) for k in layer_keys for d in dims(s[k]))
    end
    needs_domain = isempty(layer_keys) || any(d -> !(DD.name(d) in covered_names), stack_dims)
    needs_domain || return nothing

    haskey(ds, "domain") && return nothing

    attrib = Dict{String,Any}(
        "dimensions" => join((_cf_name(d) for d in stack_dims), " "),
    )
    # Aux coord names (lat/lon for ProjectedArrayLookup, inner coords for
    # MergedLookup) — same flatten logic as writevar!.
    aux_flat = String[]
    for r in stack_aux
        isnothing(r) && continue
        if r isa AbstractString
            push!(aux_flat, r)
        else
            append!(aux_flat, r)
        end
    end
    isempty(aux_flat) || (attrib["coordinates"] = join(unique(aux_flat), " "))
    grid_mapping_name = _def_grid_mapping!(ds, s)
    isnothing(grid_mapping_name) || (attrib["grid_mapping"] = grid_mapping_name)
    if any(d -> lookup(d) isa GeometryLookup, stack_dims)
        attrib["geometry"] = "geometry_container"
    end
    CDM.defVar(ds, "domain", fill('\0'), (); attrib)
    return nothing
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
    # _def_dim_var! returns aux coord var names per dim. Each writer returns
    # either nothing (no aux coords), a single name (ProjectedArrayLookup),
    # or a vector (MergedLookup, which can have several inner coords).
    # Flatten to a single list that goes in the layer's `coordinates` attr.
    aux_coord_names = String[]
    for r in _def_dim_var!(ds, A)
        isnothing(r) && continue
        if r isa AbstractString
            push!(aux_coord_names, r)
        else
            append!(aux_coord_names, r)
        end
    end
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

    # Don't write an empty-string _FillValue: it's the auto-picked sentinel
    # for String layers (see `_type_missingval1(::Type{<:AbstractString})`)
    # and would mask all legitimate empty strings as `missing` on read. The
    # inner sentinel stays "" so any actual `missing` values still encode,
    # but readers won't treat real "" as missing.
    if !isnothing(missingval_pair[1]) && !(missingval_pair[1] isa AbstractString && isempty(missingval_pair[1]))
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
    if !isempty(aux_coord_names)
        attrib["coordinates"] = join(unique(aux_coord_names), " ")
    end
    # Add grid_mapping for CRS
    grid_mapping_name = _def_grid_mapping!(ds, A)
    if !isnothing(grid_mapping_name)
        attrib["grid_mapping"] = grid_mapping_name
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

function _def_dim_var!(ds::AbstractDataset, A)
    alldims = dims(A)
    # First pass: ensure every NetCDF dimension exists, so any lookup writer
    # can reference any other dim by name (e.g. 2D aux coord vars).
    for d in alldims
        dimname = _cf_name(d)
        haskey(ds.dim, dimname) || CDM.defDim(ds, dimname, length(d))
    end
    # Second pass: hand off to the lookup-specific writer. Writers must be
    # idempotent (guard each defVar call with haskey) because shared dims in
    # a stack are visited once per layer.
    return map(d -> _def_lookup_var!(ds, d, _cf_name(d), alldims), alldims)
end
_def_lookup_var!(ds::AbstractDataset, dim::Dimension{<:AbstractNoLookup}, dimname, alldims=()) = nothing
# Generic fallback for unhandled array-shaped lookups
_def_lookup_var!(ds::AbstractDataset, dim::Dimension{<:LA.AbstractArrayLookup}, dimname, alldims=()) = nothing

# ProjectedArrayLookup: a 1D outer dim with a 2D coordinate matrix shared with
# the paired spatial dim (CF auxiliary coordinate). Writes:
#   * the 1D dim coord var (if `data` carries real values rather than just a
#     placeholder integer range)
#   * a 2D aux coord var (`lon` for X dims, `lat` for Y dims) shaped by the
#     paired (X, Y) NetCDF dimensions
# Returns the aux coord var name so writevar! can list it in the layer's
# `coordinates` attribute (which is what triggers the unaligned detection on
# load).
function _def_lookup_var!(ds::AbstractDataset, dim::Dimension{<:ProjectedArrayLookup}, dimname, alldims=())
    l = lookup(dim)
    # 1D dim coord var only when `data` is real coordinate values, and only
    # if we haven't already written it (idempotent across stack layers).
    if !(l.data isa Base.OneTo) && !haskey(ds, dimname)
        attrib = _attribdict(metadata(l))
        _cdm_set_axis_attrib!(attrib, dim)
        CDM.defVar(ds, dimname, collect(l.data), (dimname,); attrib)
    end
    # 2D aux coord vars sharing both spatial dims. The lookup carries both
    # matrices (primary in `matrix`, the other in `paired_matrix`) along with
    # the CF identity of the primary (`aux_name`). We write both lon and lat
    # vars here so the on-disk file has the full pair, regardless of whether
    # the outer dim is X/Y or plain `Dim` (CF 7.2).
    x_dim = DD.dims(alldims, XDim)
    y_dim = DD.dims(alldims, YDim)
    if isnothing(x_dim) || isnothing(y_dim)
        # Both outer dims are plain `Dim` (CF 7.2). Use them in declaration
        # order for the matrix axes - the matrix shape determines which
        # dim is which.
        plain_dims = filter(d -> !(d isa XDim || d isa YDim), alldims)
        length(plain_dims) >= 2 || return nothing
        x_dimname, y_dimname = _cf_name(plain_dims[1]), _cf_name(plain_dims[2])
    else
        x_dimname, y_dimname = _cf_name(x_dim), _cf_name(y_dim)
    end
    aux_var_name = _aux_role(l) === :lon ? "lon" : "lat"
    paired_var_name = aux_var_name == "lon" ? "lat" : "lon"
    # Write the primary
    if !haskey(ds, aux_var_name)
        # `l.dims` is set to (Lat, Lon) by the loader, but the NearestNeighbors
        # extension's `format_unaligned` rebuilds it as basedims(dims) -
        # (X{Colon}, Y{Colon}) with no lookup attached - so metadata may be
        # unavailable here. Fall back to an empty attribute dict in that case.
        primary_md_dim = aux_var_name == "lon" ? l.dims[2] : l.dims[1]
        attr = _aux_attrib(primary_md_dim)
        get!(attr, "standard_name", aux_var_name == "lon" ? "longitude" : "latitude")
        get!(attr, "units", aux_var_name == "lon" ? "degrees_east" : "degrees_north")
        CDM.defVar(ds, aux_var_name, l.matrix, (x_dimname, y_dimname); attrib=attr)
    end
    # Write the paired matrix if the lookup carries it and it isn't on disk.
    if !isnothing(l.paired_matrix) && !haskey(ds, paired_var_name)
        paired_md_dim = paired_var_name == "lon" ? l.dims[2] : l.dims[1]
        attr = _aux_attrib(paired_md_dim)
        get!(attr, "standard_name", paired_var_name == "lon" ? "longitude" : "latitude")
        get!(attr, "units", paired_var_name == "lon" ? "degrees_east" : "degrees_north")
        CDM.defVar(ds, paired_var_name, l.paired_matrix, (x_dimname, y_dimname); attrib=attr)
    end
    return isnothing(l.paired_matrix) ? aux_var_name : [aux_var_name, paired_var_name]
end

# CF aux-coord role for a ProjectedArrayLookup. Reads the `aux_name` field
# the loader set from CF metadata. Falls back to X/Y dispatch on the inner
# `dim` only when the field is `:auto` (e.g. a user-constructed lookup).
_aux_role(l::ProjectedArrayLookup) = l.aux_name === :auto ? _aux_role_from_dim(l.dim) : l.aux_name
_aux_role_from_dim(::XDim) = :lon
_aux_role_from_dim(::Any) = :lat

# Pull attrib dict from an inner aux-coord dim. Inner dims may have been
# stripped to a bare Colon-valued dim by `format_unaligned`, in which case
# `metadata` would fail - return an empty dict so the caller can fill in
# the required CF standard_name/units.
_aux_attrib(d::Dimension{<:Lookup}) = _attribdict(metadata(d))
_aux_attrib(d::Dimension) = Dict{String,Any}()

# MergedLookup: a single outer dim that bundles multiple inner coordinate
# variables sharing that one dimension (CF examples 5.3, 6.1.2, 6.2, 7.3, 7.4).
# Writes one var per inner dim, all with the merged dim as their NetCDF
# dimension:
#   * If an inner dim's name matches the merged dim's name, it becomes the
#     dim coord var (CF 6.2 sigma case).
#   * Otherwise it becomes an auxiliary coordinate var that the layer must
#     reference via the `coordinates` attribute - returned to writevar!.
function _def_lookup_var!(ds::AbstractDataset, dim::Dimension{<:MergedLookup}, dimname, alldims=())
    l = lookup(dim)
    aux_names = String[]
    for inner in l.dims
        innername = _cf_name(inner)
        if innername != dimname
            push!(aux_names, innername)
        end
        haskey(ds, innername) && continue
        inner_lookup = lookup(inner)
        attrib = _attribdict(metadata(inner_lookup))
        _maybe_write_bounds!(ds, attrib, inner, innername, dimname)
        CDM.defVar(ds, innername, collect(inner_lookup), (dimname,); attrib)
    end
    return aux_names
end

# CF section 7.1: any coord var with sampling=Intervals carries a `bounds`
# attribute pointing at a (2, N) bounds matrix. Compute the bounds from the
# in-memory lookup and write the matrix under the name the loader expects
# (preserved in metadata as `bounds`). Used by both the standalone Sampled
# writer (via the same pattern below) and the MergedLookup inner-dim writer
# so neither leaves a dangling `bounds = "..."` attribute on disk.
function _maybe_write_bounds!(ds, attrib, dim, dimname, vardim)
    l = lookup(dim)
    sampling(l) isa Intervals || return nothing
    md = metadata(l)
    # `Metadata` supports `haskey/[]` but is not `<:AbstractDict`, and the
    # source dict has String keys. The previous code used `get(md, :bounds, ...)`
    # which always missed the source name and silently wrote under e.g.
    # "a_bnds" instead of "A_bnds".
    boundskey = if md isa Lookups.Metadata && haskey(md, "bounds")
        md["bounds"]
    else
        string(dimname, "_bnds")
    end
    haskey(ds, boundskey) && (attrib["bounds"] = boundskey; return nothing)
    bounds = Dimensions.dim2boundsmatrix(dim)
    attrib["bounds"] = boundskey
    CDM.defVar(ds, boundskey, bounds, ("bnds", vardim))
    return nothing
end
function _def_lookup_var!(ds::AbstractDataset, dim::Dimension{<:GeometryLookup}, dimname, alldims=())
    # Idempotent: geometry_container is written below; if it already exists
    # the rest of the geometry vars do too.
    haskey(ds, "geometry_container") && return nothing
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
# Cyclic: written as a normal dim coord var with a `climatology` attribute
# pointing at the bounds matrix. Reverses what `_climatology_dim` did on load:
# it subtracted `duration` from each upper bound, so we add it back to recover
# the original on-disk matrix. `duration` is recoverable as
# bounds(l)[2] - val(span(l))[2, end] - the loader captured the un-shifted
# bounds tuple before mutating the matrix, exactly so this round-trip works.
function _def_lookup_var!(ds::AbstractDataset, dim::Dimension{<:Lookups.AbstractCyclic}, dimname, alldims=())
    haskey(ds, dimname) && return nothing
    l = lookup(dim)
    attrib = _attribdict(metadata(l))
    _cdm_set_axis_attrib!(attrib, dim)

    sp = span(l)
    if sp isa Explicit
        bnds_in_memory = val(sp)
        original_bounds = Lookups.bounds(l)
        duration = original_bounds[2] - bnds_in_memory[2, end]
        on_disk_bounds = copy(bnds_in_memory)
        on_disk_bounds[2, :] .+= duration
        boundskey = get(metadata(l), :climatology, "climatology_bounds")
        attrib["climatology"] = boundskey
        CDM.defVar(ds, boundskey, on_disk_bounds, ("bnds", dimname))
    end
    CDM.defVar(ds, dimname, collect(l), (dimname,); attrib)
    return nothing
end
function _def_lookup_var!(ds::AbstractDataset, dim::Dimension{<:Union{Sampled,AbstractCategorical}}, dimname, alldims=())
    # Idempotent across stack layers that share this dim
    haskey(ds, dimname) && return nothing
    # Shift lookup before conversion to Mapped
    dim = _cdm_shiftlocus(dim)
    if dim isa Y || dim isa X
        dim = convertlookup(Mapped, dim)
    end
    # Attributes
    attrib = _attribdict(metadata(dim))
    _cdm_set_axis_attrib!(attrib, dim)
    _maybe_write_bounds!(ds, attrib, dim, dimname, dimname)
    CDM.defVar(ds, dimname, collect(lookup(dim)), (dimname,); attrib)

    return nothing
end
function _def_lookup_var!(ds::AbstractDataset, dim::Dimension{<:AbstractProjected}, dimname, alldims=())
    haskey(ds, dimname) && return nothing
    # Projected lookups have CRS info - write as coordinate variable
    # The CRS itself is handled by _def_grid_mapping!
    dim = _cdm_shiftlocus(dim)
    # Attributes
    attrib = _attribdict(metadata(dim))
    _cdm_set_axis_attrib!(attrib, dim)
    _maybe_write_bounds!(ds, attrib, dim, dimname, dimname)
    CDM.defVar(ds, dimname, collect(lookup(dim)), (dimname,); attrib)
    return nothing
end

# Write one grid_mapping variable under `varname`. Returns the name or
# `nothing` if the CRS can't be serialised. Idempotent across stack layers.
function _def_grid_mapping_var!(ds::AbstractDataset, c, varname::AbstractString)
    haskey(ds, varname) && return varname
    cf = try
        convert(CFProjection, c)
    catch
        return _def_grid_mapping_wkt!(ds, c, varname)
    end
    attrib = Dict{String,Any}(cf.params)
    CDM.defVar(ds, varname, Int8(0), (); attrib)
    return varname
end

function _def_grid_mapping_wkt!(ds::AbstractDataset, c, varname::AbstractString="crs")
    haskey(ds, varname) && return varname
    wkt = try
        GeoFormatTypes.val(convert(WellKnownText2, c))
    catch
        return nothing
    end
    attrib = Dict{String,Any}(
        "crs_wkt" => wkt,
        "grid_mapping_name" => "latitude_longitude"
    )
    CDM.defVar(ds, varname, Int8(0), (); attrib)
    return varname
end

# Top-level `grid_mapping` attribute string for a layer or stack. When the
# spatial dims carry a `mappedcrs` distinct from `crs` (CF 5.10), emits the
# extended form `"crs: x y crsmapped: lat lon"` referencing two grid_mapping
# variables - one for the primary projection scoped to the dim coords, one
# for the secondary scoped to the aux coords. Otherwise emits a single name.
function _def_grid_mapping!(ds::AbstractDataset, A)
    c = crs(A)
    isnothing(c) && return nothing

    # Collect mappedcrs across spatial dims. They should all agree (it's the
    # CRS of the shared aux coords).
    mapped = nothing
    primary_coords = String[]
    mapped_coords = String[]
    for d in dims(A)
        l = lookup(d)
        l isa AbstractNoLookup && continue
        isnothing(crs(l)) && continue
        push!(primary_coords, _cf_name(d))
        if l isa ProjectedArrayLookup
            mc = l.mappedcrs
            if !isnothing(mc) && mc != crs(l)
                mapped = mc
                push!(mapped_coords, l.aux_name == :lon ? "lon" : "lat")
                isnothing(l.paired_matrix) || push!(mapped_coords, l.aux_name == :lon ? "lat" : "lon")
            end
        end
    end

    if isnothing(mapped) || isempty(mapped_coords)
        # Simple single-mapping case
        return _def_grid_mapping_var!(ds, c, "crs")
    end

    # Extended `grid_mapping = "crs: x y crsmapped: lat lon"` form.
    primary_name = _def_grid_mapping_var!(ds, c, "crs")
    mapped_name = _def_grid_mapping_var!(ds, mapped, "crsmapped")
    (isnothing(primary_name) || isnothing(mapped_name)) && return primary_name
    return string(primary_name, ": ", join(unique(primary_coords), " "), " ",
                  mapped_name, ": ", join(unique(mapped_coords), " "))
end

_cf_crs(ds, grid_mapping_key::Nothing) = nothing
function _cf_crs(ds, grid_mapping_key::AbstractString)
    if occursin(' ', grid_mapping_key)
        crss = _split_cf_attribute(grid_mapping_key)
        return map(crss) do (crs, coords)
            _cf_crs(ds, crs) => coords
        end
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
    elseif haskey(grid_mapping_attrib, "grid_mapping_name")
        # Use CFCoordinateReferenceSystems for all grid_mapping_name types
        cf = CFProjection(Dict{String,Any}(grid_mapping_attrib))
        try
            convert(WellKnownText2, cf)
        catch e
            @warn "Failed to convert CF grid mapping to CRS" exception=e grid_mapping_name=grid_mapping_attrib["grid_mapping_name"]
            nothing
        end
    else
        nothing
    end
end

function _split_cf_attribute(attribute::String)
    items = _split_attribute(attribute)
    starts = findall(contains(':'), items)
    stops = append!(starts[2:end] .- 1, length(items))
    map(starts, stops) do s, p
        items[s][1:end-1] => items[s+1:p]
    end
end

function _grid_mapping_keys(grid_mapping_key::String)
    if occursin(' ', grid_mapping_key)
        first.(_split_cf_attribute(grid_mapping_key))
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
