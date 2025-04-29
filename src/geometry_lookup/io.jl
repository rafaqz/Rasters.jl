


function _geometry_cf_encode(::GeometryOpsCore.TraitTarget{<: Union{<: GI.LineStringTrait, <: GI.MultiLineStringTrait}}, geoms)
    error("Not implemented yet")
end

function _geometry_cf_encode(::GeometryOpsCore.TraitTarget{<: Union{<: GI.PointTrait, <: GI.MultiPointTrait}}, geoms)
    error("Not implemented yet")
end



function _geometry_cf_encode(::GeometryOpsCore.TraitTarget{<: Union{<: GI.PolygonTrait, <: GI.MultiPolygonTrait}}, geoms)

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
        node_count_vec[i] = this_geom_npoints

        # push individual components of the ring
        for poly in GO.flatten(GI.PolygonTrait, geom)
            exterior_ring = GI.getexterior(poly)
            for point_idx in 1:GI.npoint(exterior_ring)-1
                point = GI.getpoint(exterior_ring, point_idx)
                xs[current_xy_index] = GI.x(point)
                ys[current_xy_index] = GI.y(point)
                current_xy_index += 1
            end
            part_node_count_vec[current_ring_index] = GI.npoint(exterior_ring)
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

    return xs, ys, node_count_vec, part_node_count_vec, interior_ring_vec

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

function _geometry_cf_decode(::Union{GI.PolygonTrait, GI.MultiPolygonTrait}, xs, ys, node_count_vec, part_node_count_vec, interior_ring_vec)
    current_xy_index = 1
    current_ring_index = 1
    current_geom_index = 1

    geoms = Vector{GO.GeometryBasics.MultiPolygon{2, Float64}}(undef, length(node_count_vec))

    for (i, npoints) in enumerate(node_count_vec)

        this_geom_npoints = npoints
        # this_geom_nholes = 
    end
end

