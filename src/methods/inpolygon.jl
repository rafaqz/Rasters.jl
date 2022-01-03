"""
    inpolygon(points, poly)

Check if a point or iterable of points is inside a polygon.

This algorithm is very efficient for many points, less so a single point.

# Arguments

- `points`: an `AbstractVector` or a `Tuple` or `Real`, Or a `Vector` of these.
- `poly`: an `AbstractVector` or nested `AbstractVector` with an inner
    `AbstractVector` or `Tuple` of `Real`. It can also be a `GeoInterface.AbstractGeometry`.

Returns a `Bool` or `AbstractVector{Bool}`.
"""
function inpolygon end
function inpolygon(points, geom::GI.AbstractGeometry; kw...)
    throw(ArgumentError("`inpolygon` works with GeoInterface polygons, but not for $(typeof(geom))"))
end
function inpolygon(points, poly::Union{GI.AbstractPolygon,GI.AbstractMultiPolygon}; kw...)
    inpolygon(points, GI.coordinates(poly); kw...)
end
function inpolygon(point::Union{Tuple,AbstractVector{<:Real}}, poly::AbstractVector; kw...)
    inpolygon([point], poly; kw...)
end
function inpolygon(points::GI.AbstractGeometry, poly::AbstractVector; kw...)
    inpolygon(GI.coordinates(points), poly, kw...)
end
function inpolygon(points::AbstractVector, poly::AbstractVector; kw...)
    inpolygon(collect(_flat_nodes(points)), poly; kw...)
end
function inpolygon(points::Union{Base.Generator,AbstractVector{<:Union{<:Tuple,<:AbstractVector{<:Real}}}}, poly::AbstractVector; kw...)
    edges = to_edges(poly)
    nodes = collect(_flat_nodes(poly))
    inpoly_matrix = inpoly2(points, nodes, edges; kw...)
    return view(inpoly_matrix, :, 1)
end
