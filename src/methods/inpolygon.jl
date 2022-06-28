"""
    inpolygon(points, poly)

Check if a point or points are inside a polygon.

This algorithm is very efficient for many points, less so a single point.

# Arguments

- `points`: an `AbstractGeom` `AbstractVector` 
- `poly`: an `AbstractVector` or nested `AbstractVector` with an inner
    `AbstractVector` or `Tuple` of `Real`, or a `GeoInterface.AbstractPolygon`.

Returns a `Bool` or `AbstractVector{Bool}`.
"""
inpolygon(points, geom; kw...) = _inpolygon(points, geom; kw...)

_inpolygon(points, geom; kw...) = _inpolygon(GI.trait(points), points, geom; kw...)
function _inpolygon(::Nothing, points::AbstractVector, geom; kw...)
    edges, nodes = to_edges_and_nodes(geom)
    inpoly_matrix = inpoly2(points, nodes, edges; kw...)
    return view(inpoly_matrix, :, 1)
end
function _inpolygon(::Nothing, points, geom; kw...)
    throw(ArgumentError("object is not an GeoInterface geometry or feature."))
end
function _inpolygon(::GI.AbstractFeatureTrait, points, geom; kw...)
    _inpolygon(GI.geometry(points), geom; kw...)
end
function _inpolygon(::GI.AbstractFeatureCollectionTrait, points, geom; kw...)
    [_inpolygon(f, geom, kw...) for f in GI.getfeature(points)]
end
function _inpolygon(pt::GI.AbstractGeometryTrait, points, geom; kw...)
    _inpolygon(pt, points, GI.trait(geom), geom; kw...)
end
function _inpolygon(pt::GI.AbstractGeometryTrait, points, gt::GI.AbstractFeatureTrait, feature; kw...)
    _inpolygon(pt, points, GI.geometry(feature); kw...)
end
function _inpolygon(pt::GI.AbstractGeometryTrait, points, gt::GI.AbstractFeatureCollectionTrait, fc; kw...)
    matches = map(f -> _inpolygon(pt, points, f; kw...), GI.getfeature(fc))
    any.(matches)
end
function _inpolygon(::GI.AbstractGeometryTrait, points, ::GI.AbstractGeometryTrait, geom; kw...)
    edges, nodes = to_edges_and_nodes(geom)
    tuplepoints = [(GI.x(p), GI.y(p)) for p in GI.getpoint(points)]
    inpoly_matrix = inpoly2(tuplepoints, nodes, edges; kw...)
    return all(view(inpoly_matrix, :, 1))
end
function _inpolygon(::GI.AbstractPointTrait, points, ::GI.AbstractGeometryTrait, geom; kw...)
    edges, nodes = to_edges_and_nodes(geom)
    return first(inpoly2([(GI.x(points), GI.y(points))], nodes, edges; kw...))
end

# Copied from PolygonInbounds, to add extra keyword arguments
# PR to include these when this has solidified
function inpoly2(vert, node, edge=zeros(Int);
    atol::T=0.0, rtol::T=NaN, iyperm=nothing,
    vmin=nothing, vmax=nothing, pmin=nothing, pmax=nothing
) where T<:AbstractFloat
    rtol = !isnan(rtol) ? rtol : iszero(atol) ? eps(T)^0.85 : zero(T)
    poly = PolygonInbounds.PolygonMesh(node, edge)
    points = PolygonInbounds.PointsInbound(vert)

    vmin = isnothing(vmin) ? minimum(points) : vmin
    vmax = isnothing(vmax) ? maximum(points) : vmax
    pmin = isnothing(pmin) ? minimum(poly) : pmin
    pmax = isnothing(pmax) ? maximum(poly) : pmax

    lbar = sum(pmax - pmin)
    tol = max(abs(rtol * lbar), abs(atol))

    ac = PolygonInbounds.areacount(poly)
    stat = ac > 1 ? falses(length(points), 2, ac) : falses(length(points), 2)
    # flip coordinates so expected effort is minimal
    dvert = vmax .- vmin
    if isnothing(iyperm)
        ix = dvert[1] < dvert[2] ? 1 : 2
        iyperm = sortperm(points, 3 - ix)
    else
        ix = 1
    end

    PolygonInbounds.inpoly2!(points, iyperm, poly, ix, tol, stat)
    return stat
end
