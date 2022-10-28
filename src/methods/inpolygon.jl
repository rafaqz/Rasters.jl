"""
    inpolygon(points, poly)

Check if a point or points are inside a polygon.

This algorithm is efficient for many points, less so for a single point.

# Arguments

- `points`: a GeoInterface.jl compatable object, or `AbstractVector` of objects.
- `poly`: a GeoInterface.jl compatable polygon, or a feature or feature collection
    of polygons.

Returns a `Bool` or `AbstractVector{Bool}`.
"""
inpolygon(points, poly; kw...) = _inpolygon(points, poly; kw...)

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


using PolygonInbounds: vertex, edgecount, flipio!, edgeindex, searchfirst

# Copied from PolygonInbounds, to add extra keyword arguments.
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
    # stat = ac > 1 ? fill(false, length(points), 1, ac) : 
    stat = fill(false, length(points), 2)
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

function inpoly2!(points, iyperm, poly, ix::Integer, veps::T, stat::S) where {N,T<:AbstractFloat,S<:AbstractArray{Bool,N}}
    nvrt = length(points)   # number of points to be checked
    nedg = edgecount(poly)  # number of edges of the polygon mesh
    vepsx = vepsy = veps
    iy = 3 - ix

    #----------------------------------- loop over polygon edges
    for epos = 1:nedg

        inod = edgeindex(poly, epos, 1)  # from
        jnod = edgeindex(poly, epos, 2)  # to
        # swap order of vertices
        if vertex(poly, inod, iy) > vertex(poly, jnod, iy)
            inod, jnod = jnod, inod
        end

        #------------------------------- calc. edge bounding-box
        xone = vertex(poly, inod, ix)
        yone = vertex(poly, inod, iy)
        xtwo = vertex(poly, jnod, ix)
        ytwo = vertex(poly, jnod, iy)

        xmin0 = min(xone, xtwo)
        xmax0 = max(xone, xtwo)
        xmin = xmin0 - vepsx
        xmax = xmax0 + vepsx
        ymin = yone - vepsy
        ymax = ytwo + vepsy

        ydel = ytwo - yone
        xdel = xtwo - xone
        xysq = xdel^2 + ydel^2
        feps = sqrt(xysq) * veps

        # find top points[:,iy] < ymin by binary search
        ilow = searchfirst(points, iy, iyperm, ymin)
        #------------------------------- calc. edge-intersection
        # loop over all points with y ∈ [ymin,ymax)
        for jpos = ilow:nvrt
            jorig = iyperm[jpos]
            ypos = vertex(points, jorig, iy)
            ypos > ymax && break 
            xpos = vertex(points, jorig, ix)

            if xpos >= xmin
                if xpos <= xmax
                    #--------- inside extended bounding box of edge
                    mul1 = ydel * (xpos - xone)
                    mul2 = xdel * (ypos - yone)
                    # if abs(mul2 - mul1) <= feps
                        #------- distance from line through edge less veps
                        # mul3 = xdel * (2xpos-xone-xtwo) + ydel * (2ypos-yone-ytwo)
                        # if abs(mul3) <= xysq ||
                        #    hypot(xpos- xone, ypos - yone) <= veps ||
                        #    hypot(xpos- xtwo, ypos - ytwo) <= veps
                        #     # ---- round boundaries around endpoints of edge
                        #     setonbounds!(poly, stat, jorig, epos)
                        # end
                            #----- left of line && ypos exact to avoid multiple counting
                    # end
                    if mul1 < mul2 && yone <= ypos < ytwo
                    elseif mul1 < mul2 && yone <= ypos < ytwo
                        #----- left of line && ypos exact to avoid multiple counting
                        flipio!(poly, stat, jorig, epos)
                    end
                end
            else # xpos < xmin - left of bounding box
                if yone <= ypos <  ytwo
                    #----- ypos exact to avoid multiple counting
                    flipio!(poly, stat, jorig, epos)
                end
            end
        end
    end
    stat
end
