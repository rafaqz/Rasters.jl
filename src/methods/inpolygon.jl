"""
    inpolygon(points, poly)

Check if a point or points are inside a polygon.

This algorithm is very efficient for many points, less so a single point.

# Arguments

- `points`: an `AbstractVector` or a `Tuple` or `Real`, an `AbstractVector` of these,
    or a `GeoInterface.AbstractGeometry`.
- `poly`: an `AbstractVector` or nested `AbstractVector` with an inner
    `AbstractVector` or `Tuple` of `Real`, or a `GeoInterface.AbstractPolygon`.

Returns a `Bool` or `AbstractVector{Bool}`.
"""
function inpolygon(point::Union{Tuple,AbstractVector{<:Real}}, poly; kw...)
    inpolygon([point], poly; kw...)[1]
end
function inpolygon(points::GI.AbstractGeometry, poly; kw...)
    inpolygon(collect(_flat_nodes(GI.coordinates(points))), poly, kw...)
end
function inpolygon(points, poly; kw...)
    vecs = _node_vecs(poly)
    edges, nodes = to_edges_and_nodes(vecs)
    inpoly_matrix = inpoly2(points, nodes, edges; kw...)
    return view(inpoly_matrix, :, 1)
end
