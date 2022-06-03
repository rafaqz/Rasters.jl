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
inpolygon(geom, poly; kw...) = inpolygon(GI.geomtrait(points), geom, poly; kw...) 
function inpolygon(::Nothing, geom, poly; kw...)
    edges, nodes = to_edges_and_nodes(poly)
    inpoly_matrix = inpoly2(points, nodes, edges; kw...)
    return view(inpoly_matrix, :, 1)
end
function inpolygon(::GI.AbstractFeatureTrait, feature, poly; kw...)
    geom = GI.geometry(feature)
    inpolygon(GI.geomtrait(geom), geom, kw...)
end
function inpolygon(t::GI.AbstractGeometryTrait, geom, poly; kw...)
    return inpolygon(GI.getpoint(t, geom), poly, kw...)
end
function inpolygon(::GI.AbstractPointTrait, geom, poly; kw...)
    inpolygon((GI.x(point), GI.y(point)), poly; kw...)[1]
end
