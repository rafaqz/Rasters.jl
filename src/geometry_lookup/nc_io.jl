using NCDatasets
import NCDatasets.CommonDataModel as CDM

import Rasters


ds = NCDataset("/Users/anshul/i.nc")
var = ds["someData"]

knowndims = Rasters._dims(var)  

unknowndims_idxs = findall(Rasters.isnolookup ∘ Rasters.lookup, knowndims)

if length(unknowndims_idxs) > 1
    throw(ArgumentError("Only one unknown dimension is supported"))
elseif length(unknowndims_idxs) == 0
    return knowndims
end

u_idx = only(unknowndims_idxs)

u_dim_name = CDM.dimnames(var)[u_idx]

has_geometry = haskey(CDM.attribs(var), "geometry")
if !has_geometry
    throw(ArgumentError("Variable $u_dim_name does not have a geometry attribute"))
end
geometry_container_varname = CDM.attribs(var)["geometry"]
geometry_container_var = ds[geometry_container_varname]

geometry_container_attribs = CDM.attribs(geometry_container_var)

haskey(geometry_container_attribs, "geometry_type") && 
geometry_container_attribs["geometry_type"] == "polygon" ||
throw(ArgumentError("We only support polygon geometry types at this time, got $geometry_type"))

@assert haskey(ds, geometry_container_attribs["node_count"])
node_count_var = ds[geometry_container_attribs["node_count"]]
only(CDM.dimnames(node_count_var)) != u_dim_name && throw(ArgumentError("node_count variable $u_dim_name does not match the unknown dimension $u_dim_name"))

node_count = collect(node_count_var)
node_coordinates = collect(zip(getindex.((ds,), split(geometry_container_attribs["node_coordinates"], " "))...))
part_node_count = collect(ds[geometry_container_attribs["part_node_count"]])
interior_ring = collect(ds[geometry_container_attribs["interior_ring"]])

current_xy_index = 1
current_ring_index = 1
current_geom_index = 1


# Initialize variables for ring assembly
start = 1
stop = part_node_count[1]
rings = [node_coordinates[start:stop]]

# Assemble all rings
geoms = Vector{GI.MultiPolygon{2, Float64}}(undef, length(node_count))

for i in 2:length(part_node_count)
    start = stop + 1
    stop = start + part_node_count[i] - 1
    push!(rings, node_coordinates[start:stop])
    # Ensure rings are closed by adding the first point at the end
    push!(rings[end], rings[end][1])
end

# Assemble multipolygons
current_ring = 1
for (geom_idx, total_nodes) in enumerate(node_count)
    # Find all rings that belong to this polygon
    polygon_rings = []
    accumulated_nodes = 0
    
    while current_ring <= length(part_node_count) && accumulated_nodes < total_nodes
        ring = rings[current_ring]
        push!(polygon_rings, (GI.LinearRing(ring), interior_ring[current_ring]))
        accumulated_nodes += part_node_count[current_ring]
        current_ring += 1
    end
    
    # Create polygons from rings
    polygons = GI.Polygon[]
    current_exterior = nothing
    current_holes = GI.LinearRing[]
    
    for (ring, is_interior) in polygon_rings
        if is_interior == 0
            # If we have a previous exterior ring, create a polygon with it and its holes
            if !isnothing(current_exterior)
                push!(polygons, GI.Polygon([current_exterior, current_holes...]))
                current_holes = GI.LinearRing[]
            end
            current_exterior = ring
        else
            push!(current_holes, ring)
        end
    end
    
    # Add the last polygon if we have an exterior ring
    if !isnothing(current_exterior)
        push!(polygons, GI.Polygon([current_exterior, current_holes...]))
    end
    
    # Create multipolygon from all polygons
    geoms[geom_idx] = GI.MultiPolygon(polygons)
end

