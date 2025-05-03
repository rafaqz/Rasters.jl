#=
# Geometry encoding and decoding from CF conventions

Encode functions will always return a named tuple with the standard 
=#
function _geometry_cf_encode(::Union{GI.PointTrait, GI.MultiPointTrait}, geoms)
    if all(x -> GI.trait(x) isa GI.PointTrait, geoms)
        return (;
            node_coordinates_x = GI.x.(geoms),
            node_coordinates_y = GI.y.(geoms),
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

    return (;
        node_coordinates_x = flat_xs,
        node_coordinates_y = flat_ys,
        node_count = GI.npoint.(geoms)
    )
end


function _geometry_cf_encode(::Union{GI.LineStringTrait, GI.MultiLineStringTrait}, geoms)
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

    # If all geoms are linestrings, we can take a fast path.
    if all(x -> GI.trait(x) isa GI.LineStringTrait, geoms)
        return (;
            node_coordinates_x = flat_xs,
            node_coordinates_y = flat_ys,
            node_count = GI.npoint.(geoms)
        )
    end

    # Otherwise, we need to encode part_node_count for multilinestrings.
    return (;
        node_coordinates_x = flat_xs,
        node_coordinates_y = flat_ys,
        part_node_count = collect(GO.flatten(GI.npoint, GI.LineStringTrait, geoms)),
        node_count = GI.npoint.(geoms)
    )
end

function _geometry_cf_encode(::Union{GI.PolygonTrait, GI.MultiPolygonTrait}, geoms)

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
            # Main.@infiltrate
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
    return (;
        node_coordinates_x = xs, 
        node_coordinates_y = ys, 
        node_count = node_count_vec, 
        part_node_count = part_node_count_vec, 
        interior_ring = interior_ring_vec
    )
end

#=
function _def_dim_var!(ds::AbstractDataset, dim::Dimension{<: GeometryLookup})
    dimname = lowercase(string(DD.name(dim)))
    haskey(ds.dim, dimname) && return nothing
    CDM.defDim(ds, dimname, length(dim))
    lookup(dim) isa NoLookup && return nothing
    attribdict = _attribdict(NoMetadata())


    CDM.defVar(ds, dimname, Vector(index(dim)), (dimname,); attrib=attrib)
    
end
=#


function _geometry_cf_decode(::Union{GI.PolygonTrait, GI.MultiPolygonTrait}, ds, geometry_container_attribs; crs = nothing)
    # First of all, we assert certain things about the geometry container and what it has.
    @assert haskey(ds, geometry_container_attribs["node_count"])
    node_count_var = ds[geometry_container_attribs["node_count"]]
    # only(CDM.dimnames(node_count_var)) != u_dim_name && throw(ArgumentError("node_count variable $u_dim_name does not match the unknown dimension $u_dim_name"))

    # Load and create all the data we need.
    node_count = collect(node_count_var)
    node_coordinates = collect(zip(getindex.((ds,), split(geometry_container_attribs["node_coordinates"], " "))...))

    # We can take a fast path for polygons, if we know that there are no multipart polygons.
    if !haskey(geometry_container_attribs, "part_node_count")
        node_count_stops = cumsum(node_count)
        node_count_starts = [1, node_count_stops[1:end-1] .+ 1...]
        return map(node_count_starts, node_count_stops) do start, stop
            GI.Polygon([GI.LinearRing(node_coordinates[start:stop]; crs)]; crs)
        end
    end


    part_node_count = collect(ds[geometry_container_attribs["part_node_count"]])
    interior_ring = collect(ds[geometry_container_attribs["interior_ring"]])

    # First, we assemble all the rings.  That's the slightly complex part.
    # After rings are assembled, we assemble the polygons and multipolygons from the rings.

    # Initialize variables for ring assembly
    start = 1
    stop = part_node_count[1]
    rings = [node_coordinates[start:stop]]
    push!(rings[end], node_coordinates[start])

    # Assemble all rings
    for i in 2:length(part_node_count)
        start = stop + 1
        stop = start + part_node_count[i] - 1
        push!(rings, node_coordinates[start:stop])
        # Ensure rings are closed by adding the first point at the end
        push!(rings[end], node_coordinates[start])
    end

    # Now, we proceed to assemble the polygons and multipolygons from the rings.
    # TODO: no better way to get the tuple type, at least for now.
    _lr = GI.LinearRing([(0.0, 0.0), (1.0, 1.0), (1.0, 0.0), (0.0, 0.0)]; crs)
    _p = GI.Polygon([_lr]; crs)
    _mp = GI.MultiPolygon([_p]; crs)
    geoms = Vector{typeof(_mp)}(undef, length(node_count))
    # Assemble multipolygons
    current_ring = 1
    for (geom_idx, total_nodes) in enumerate(node_count)
        # Find all rings that belong to this polygon
        polygon_rings = Tuple{typeof(_lr), Int8}[]
        accumulated_nodes = 0
        
        while current_ring <= length(part_node_count) && accumulated_nodes < total_nodes
            ring = rings[current_ring]
            push!(polygon_rings, (GI.LinearRing(ring; crs), interior_ring[current_ring]))
            accumulated_nodes += part_node_count[current_ring]
            current_ring += 1
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



function _geometry_cf_decode(::Union{GI.LineStringTrait, GI.MultiLineStringTrait}, ds, geometry_container_attribs; crs = nothing)
    @assert haskey(ds, geometry_container_attribs["node_count"])
    node_count_var = ds[geometry_container_attribs["node_count"]]

    # Load and create all the data we need.
    node_count = collect(node_count_var)
    node_coordinates = collect(zip(getindex.((ds,), split(geometry_container_attribs["node_coordinates"], " "))...))

    # we can use a fast path for lines, if we know that there are no multipart lines.
    if !haskey(geometry_container_attribs, "part_node_count")
        node_count_stops = cumsum(node_count)
        node_count_starts = [1, node_count_stops[1:end-1] .+ 1...]
        return GI.LineString.(getindex.((node_coordinates,), (:).(node_count_starts, node_count_stops)); crs)
    end

    # otherwise, we need to decode the multipart lines.
    part_node_count = collect(ds[geometry_container_attribs["part_node_count"]])

    # Initialize variables for line assembly
    start = 1
    stop = part_node_count[1]
    lines = [node_coordinates[start:stop]]

    # Assemble all lines
    for i in 2:length(part_node_count)
        start = stop + 1
        stop = start + part_node_count[i] - 1
        push!(lines, node_coordinates[start:stop])
    end

    # Now assemble the multilinestrings
    _ls = GI.LineString(node_coordinates[1:2]; crs)
    _mls = GI.MultiLineString([_ls]; crs)
    geoms = Vector{typeof(_mls)}(undef, length(node_count))

    # Assemble multilinestrings
    current_line = 1
    for (geom_idx, total_nodes) in enumerate(node_count)
        # Find all lines that belong to this multilinestring
        multilinestring_lines = typeof(_ls)[]
        accumulated_nodes = 0
        
        while current_line <= length(part_node_count) && accumulated_nodes < total_nodes
            line = lines[current_line]
            push!(multilinestring_lines, GI.LineString(line; crs))
            accumulated_nodes += part_node_count[current_line]
            current_line += 1
        end
        
        # Create multilinestring from all lines
        geoms[geom_idx] = GI.MultiLineString(multilinestring_lines; crs)
    end

    return geoms
end

function _geometry_cf_decode(::Union{GI.PointTrait, GI.MultiPointTrait}, ds, geometry_container_attribs; crs = nothing)
    
    node_coordinates = collect(zip(getindex.((ds,), split(geometry_container_attribs["node_coordinates"], " "))...))
    # We can take a fast path for points, if we know that there are no multipoints
    if haskey(geometry_container_attribs, "node_count")
        @assert haskey(ds, geometry_container_attribs["node_count"])
        node_count_var = ds[geometry_container_attribs["node_count"]]
        node_count = collect(node_count_var)
        # The code below could be a fast path, but we don't want
        # to arbitrarily change the output type of the decoder.
        # MultiPoints should always roundtrip and write as multipoints.
        # if !all(==(1), node_count)
        # do nothing
        # else
        # return a fast path 
        # end
        # we have multipoints
        node_count_stops = cumsum(node_count)
        node_count_starts = [1, node_count_stops[1:end-1] .+ 1...]
        return map(node_count_starts, node_count_stops) do start, stop
            GI.MultiPoint(node_coordinates[start:stop]; crs)
        end
    end

    # finally, if we have no node count, or all node counts are 1, we just return the points
    return GI.Point.(node_coordinates; crs)

end