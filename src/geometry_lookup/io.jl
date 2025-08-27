# Geometry encoding and decoding from CF conventions (7.5)
_cf_geometry_encode(geoms) = _cf_geometry_encode(GI.trait(first(geoms)), geoms)
_cf_geometry_encode(trait::GI.AbstractTrait, geoms) = throw(ArgumentError("Geometry trait $trait not currently handled by Rasters"))
function _cf_geometry_encode(::Union{GI.PointTrait,GI.MultiPointTrait}, geoms)
    if all(x -> GI.trait(x) isa GI.PointTrait, geoms)
        return (;
            geometry_type = "point",
            x = GI.x.(geoms),
            y = GI.y.(geoms),
        )
    end
    # else: we have some multipolygons in here
    npoints = sum(GI.npoint, geoms)
    flat_xs = Vector{Float64}(undef, npoints)
    flat_ys = Vector{Float64}(undef, npoints)

    i::Int = 0
    # flatten is guaranteed to iterate in order.
    flattener = GO.flatten(GI.PointTrait, geoms) do point
        i += 1
        flat_xs[i] = GI.x(point)
        flat_ys[i] = GI.y(point)
    end
    # iterate over flattener to populate the arrays,
    # without allocating an extra array.
    foreach(identity, flattener)

    centroids = GO.centroid.(geoms)
    return (;
        geometry_type = "line",
        x = flat_xs,
        y = flat_ys,
        lon = first.(centroids),
        lat = last.(centroids),
        node_count = GI.npoint.(geoms)
    )
end
function _cf_geometry_encode(::Union{GI.LineStringTrait, GI.MultiLineStringTrait}, geoms)
    # There is a fast path without encoding part_node_count if all geoms are linestrings.
    npoints = sum(GI.npoint, geoms)
    flat_xs = Vector{Float64}(undef, npoints)
    flat_ys = Vector{Float64}(undef, npoints)

    i::Int = 0
    # flatten is guaranteed to iterate in order.
    flattener = GO.flatten(GI.PointTrait, geoms) do point
        i += 1
        flat_xs[i] = GI.x(point)
        flat_ys[i] = GI.y(point)
    end
    # iterate over flattener to populate the arrays,
    # without allocating an extra array.
    foreach(identity, flattener)
    attribs = Dict{String,Any}("geometry_type" => "linestring")

    # If all geoms are linestrings, we can take a fast path.
    centroids = GO.centroid.(geoms)
    if all(x -> GI.trait(x) isa GI.LineStringTrait, geoms)
        return (;
            geometry_type = "line",
            x = flat_xs,
            y = flat_ys,
            lon = first.(centroids),
            lat = last.(centroids),
            node_count = GI.npoint.(geoms)
        )
    end

    # Otherwise, we need to encode part_node_count for multilinestrings.
    return (;
        geometry_type = "line",
        x = flat_xs,
        y = flat_ys,
        lon = first.(centroids),
        lat = last.(centroids),
        part_node_count = collect(GO.flatten(GI.npoint, GI.LineStringTrait, geoms)),
        node_count = GI.npoint.(geoms)
    )
end
function _cf_geometry_encode(::Union{GI.PolygonTrait, GI.MultiPolygonTrait}, geoms)
    ngeoms = length(geoms)
    nrings = GO.applyreduce(GI.nring, +, GI.PolygonTrait(), geoms; init = 0, threaded = false)
    n_points_per_geom_vec = GI.npoint.(geoms)
    total_n_points = sum(n_points_per_geom_vec) - nrings

    # Create a vector of the total number of points
    xs = fill(0.0, total_n_points)
    ys = fill(0.0, total_n_points)

    node_count_vec = fill(0, ngeoms)
    part_node_count_vec = fill(0, nrings)
    interior_ring_vec = fill(0, nrings)

    current_xy_index = 1
    current_ring_index = 1

    for (i, geom) in enumerate(geoms)
        this_geom_npoints = GI.npoint(geom)
        # Bear in mind, that the last point (which == first point) 
        # of the linear ring is removed when encoding, so not included
        # in the node count.
        node_count_vec[i] = this_geom_npoints - GI.nring(geom)

        # push individual components of the ring
        for poly in GO.flatten(GI.PolygonTrait, geom)
            exterior_ring = GI.getexterior(poly)
            for point_idx in 1:GI.npoint(exterior_ring)-1
                point = GI.getpoint(exterior_ring, point_idx)
                xs[current_xy_index] = GI.x(point)
                ys[current_xy_index] = GI.y(point)
                current_xy_index += 1
            end
            part_node_count_vec[current_ring_index] = GI.npoint(exterior_ring)-1
            interior_ring_vec[current_ring_index] = 0
            current_ring_index += 1

            if GI.nring(poly) == 1
                continue
            else
                for hole in GI.gethole(poly)
                    for point_idx in 1:GI.npoint(hole)-1
                        point = GI.getpoint(hole, point_idx)
                        xs[current_xy_index] = GI.x(point)
                        ys[current_xy_index] = GI.y(point)
                        current_xy_index += 1
                    end
                    part_node_count_vec[current_ring_index] = GI.npoint(hole)-1
                    interior_ring_vec[current_ring_index] = 1
                    current_ring_index += 1
                end
            end
        end
    end
    # Note: this does not encode the `lat` and `lon` variables that might hold a representative point of the polygon, like a centroid.
    # The names in this named tuple are standard CF conventions.
    # node_coordinates_x and node_coordinates_y are the coordinates of the nodes of the rings.
    # but cf encodes them as a space separated string.  That's the only difference.
    centroids = GO.centroid.(geoms)
    return (;
        geometry_type = "polygon",
        x = xs, 
        y = ys, 
        lon = first.(centroids),
        lat = last.(centroids),
        node_count = node_count_vec, 
        part_node_count = part_node_count_vec, 
        interior_ring = interior_ring_vec
    )
end

# CF standards Geometry decoding (7.5)
# `geometry` is the CF attributes dict from the variable linked to a 'geometry" attribute.
# `ds` is any CommonDataModel.AbstractDataset
function _cf_geometry_decode(ds::AbstractDataset, geometry; kw...)
    geometry_type = geometry["geometry_type"]
    trait = if geometry_type == "point"
        haskey(geometry, "node_count") ? GI.MultiPointTrait() : GI.PointTrait()
    elseif geometry_type == "line"
        haskey(geometry, "part_node_count") ? GI.MultiLineStringTrait() : GI.LineStringTrait()
    elseif geometry_type == "polygon"
        haskey(geometry, "part_node_count") ?  GI.MultiPolygonTrait() : GI.PolygonTrait()
    end
    return _cf_geometry_decode(trait, ds, geometry; kw...)
end
function _cf_geometry_decode(::GI.MultiPolygonTrait, ds, geometry; crs=nothing)
    rings = _split_inner_geoms(ds, geometry; autoclose=true)
    node_count = _cf_node_count(ds, geometry)
    interior_ring = _cf_interior_ring(ds, geometry)
    # Now, we proceed to assemble the polygons and multipolygons from the rings.
    # TODO: no better way to get the tuple type, at least for now.
    _lr = GI.LinearRing(first(rings); crs)
    _p = GI.Polygon([_lr]; crs)
    _mp = GI.MultiPolygon([_p]; crs)
    geoms = Vector{typeof(_mp)}(undef, length(node_count))

    # Assemble multipolygons
    current_ring = 1
    for (geom_idx, total_nodes) in enumerate(node_count)
        # Find all rings that belong to this polygon
        polygon_rings = Tuple{typeof(_lr), Int}[]
        
        n_points_added = 0
        while current_ring <= length(rings) && n_points_added < total_nodes
            ring = rings[current_ring]
            push!(polygon_rings, (GI.LinearRing(ring; crs), interior_ring[current_ring]))
            current_ring += 1
            n_points_added += length(ring)
        end
        
        # Create polygons from rings
        polygons = typeof(_p)[]
        current_exterior = nothing
        current_holes = typeof(_lr)[]
        
        for (ring, is_interior) in polygon_rings
            if is_interior == 0
                # If we have a previous exterior ring, create a polygon with it and its holes
                if !isnothing(current_exterior)
                    push!(polygons, GI.Polygon([current_exterior, current_holes...]; crs))
                    current_holes = typeof(_lr)[]
                end
                current_exterior = ring
            else
                push!(current_holes, ring)
            end
        end
        
        # Add the last polygon if we have an exterior ring
        if !isnothing(current_exterior)
            push!(polygons, GI.Polygon([current_exterior, current_holes...]; crs))
        end
        # Create multipolygon from all polygons
        geoms[geom_idx] = GI.MultiPolygon(polygons; crs)
    end
    return geoms
end
function _cf_geometry_decode(::GI.MultiLineStringTrait, ds, geometry; crs=nothing)
    node_count = _cf_node_count(ds, geometry)
    lines = _split_inner_geoms(ds, geometry; autoclose = false)

    _ls = GI.LineString(lines[1]; crs)
    _mls = GI.MultiLineString([_ls]; crs)
    geoms = Vector{typeof(_mls)}(undef, length(node_count))

    # Assemble multilinestrings
    current_line = 1
    for (geom_idx, total_nodes) in enumerate(node_count)
        # Find all lines that belong to this multilinestring
        multilinestring_lines = typeof(_ls)[]
        nodes_added = 0
        while nodes_added < total_nodes
            line = lines[current_line]
            push!(multilinestring_lines, GI.LineString(line; crs))
            current_line += 1
            nodes_added += length(line)
        end
        # Create multilinestring from all lines
        geoms[geom_idx] = GI.MultiLineString(multilinestring_lines; crs)
    end

    return geoms
end
function _cf_geometry_decode(::GI.PointTrait, ds, geometry; crs=nothing)
    node_coordinates = _cf_node_coordinates(ds, geometry)
    # Just wrap raw coordinates as Points
    return GI.Point.(node_coordinates; crs)
end
function _cf_geometry_decode(::GI.LineStringTrait, ds, geometry; crs=nothing)
    node_count = _cf_node_count(ds, geometry)
    node_coordinates = _cf_node_coordinates(ds, geometry)
    # Split coordinates to separate LineStrings by node count ranges
    return map(_node_ranges(node_count)) do range
        GI.LineString(node_coordinates[range]; crs)
    end
end
function _cf_geometry_decode(::GI.MultiPointTrait, ds, geometry; crs = nothing)
    node_count = _cf_node_count(ds, geometry)
    node_coordinates = _cf_node_coordinates(ds, geometry)
    # Split coordinates to separate MultiPoint by node count ranges
    return map(_node_ranges(node_count)) do range
        GI.MultiPoint(node_coordinates[range]; crs)
    end
end
function _cf_geometry_decode(::GI.PolygonTrait, ds, geometry; crs=nothing)
    node_count = _cf_node_count(ds, geometry)
    node_coordinates = _cf_node_coordinates(ds, geometry)
    # Split coordinates to separate single-ring Polygons by node count ranges
    return map(_node_ranges(node_count)) do range
        GI.Polygon([GI.LinearRing(node_coordinates[range]; crs)]; crs)
    end
end

function _node_ranges(node_count) 
    cum_node_count = cumsum(node_count)
    ranges = Vector{UnitRange{Int}}(undef, length(node_count))
    for i in eachindex(ranges)
        ranges[i] = i == 1 ? (1:node_count[i]) : ((cum_node_count[i-1]+1):cum_node_count[i])
    end
    return ranges
end

function _split_inner_geoms(ds, geometry; autoclose=false)
    part_node_count = _cf_part_node_count(ds, geometry)
    node_coordinates = _cf_node_coordinates(ds, geometry)
    # Initialize variables for ring assembly
    start = 1
    stop = part_node_count[1]
    rings = [node_coordinates[start:stop]]
    autoclose && push!(rings[end], node_coordinates[start])

    # Assemble all rings
    for i in 2:length(part_node_count)
        start = stop + 1
        stop = start + part_node_count[i] - 1
        push!(rings, node_coordinates[start:stop])
        # Ensure rings are closed by adding the first point at the end
        autoclose && push!(rings[end], node_coordinates[start])
    end
    return rings
end

function _cf_node_coordinates(ds, geometry) 
    coords = map(_cf_node_coordinate_names(geometry)) do coordname
        collect(CDM.variable(ds, coordname))
    end
    return collect(zip(coords...))
end
_cf_node_coordinate_names(geometry) = split(geometry["node_coordinates"], ' ')
_cf_node_count(ds, geometry) = collect(CDM.variable(ds, geometry["node_count"]))
_cf_part_node_count(ds, geometry) = collect(CDM.variable(ds, geometry["part_node_count"]))
_cf_interior_ring(ds, geometry) = collect(CDM.variable(ds, geometry["interior_ring"]))