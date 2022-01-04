const Pt{T<:Real} = Union{AbstractVector{T},NTuple{<:Any,T}}
const Poly = AbstractVector{<:Union{NTuple{<:Any,<:Real},AbstractVector{<:Real}}}

const DEFAULT_POINT_ORDER = (XDim, YDim)

# _fill_geometry!
# Fill a raster with `fill` where it interacts with a geometry.
# We reduce all geometries to vectors with `coordinates` so that a simple
# vector can also represent a polygon. In future we should do this lazily.
function _fill_geometry!(B::AbstractRaster{Bool}, coll::GI.AbstractFeatureCollection; order, kw...)
    if bbox_might_overlap(B, coll, order)
        foreach(coll.features) do feature
            _fill_geometry!(B, feature; order, kw...)
        end
    end
    return B
end
function _fill_geometry!(B::AbstractRaster{Bool}, feature::GI.AbstractFeature; order, kw...)
    if bbox_might_overlap(B, feature, order)
        _fill_geometry!(B, feature.geometry; order, kw...)
    end
    return B
end
function _fill_geometry!(B::AbstractRaster{Bool}, geom::GI.AbstractGeometry;
    shape=nothing, order=DEFAULT_ORDER, kw...
)
    order = dims(B, order)
    shape = isnothing(shape) ? _geom_shape(geom) : shape
    if bbox_might_overlap(B, geom, order)
        _fill_geometry!(B, GI.coordinates(geom); shape, order, kw...)
    end
    return B
end
function _fill_geometry!(B::AbstractRaster{Bool}, geom; shape=:polygon, order, kw...)
    geom = _flat_nodes(geom)
    gbounds = _geom_bounds(geom, order)
    abounds = bounds(dims(B, order)) # Only mask if the gemoetry bounding box overlaps the array bounding box
    bounds_overlap(gbounds, abounds) || return B
    if shape === :polygon
        _fill_polygon!(B, geom; polybounds=gbounds, order, kw...)
    elseif shape === :line
        _fill_linestring!(B, geom; order, kw...)
    elseif shape === :point
        _fill_point!(B, geom; order, kw...)
    else
        _shape_error(shape)
    end
    return B
end

# _fill_polygon!
# Fill a raster with `fill` where pixels are inside a polygon
# `boundary` determines how edges are handled
function _fill_polygon!(B::AbstractRaster{Bool}, poly::Base.Iterators.Flatten; kw...)
    _fill_polygon!(B, collect(poly); kw...)
end
function _fill_polygon!(B::AbstractRaster{Bool}, poly; polybounds, order, fill=true, boundary=:center, kw...)
    # TODO take a view of B for the polygon
    # We need a tuple of all the dims in `order`
    # We also need the index locus to be the center so we are
    # only selecting cells more than half inside the polygon
    shifted_dims = map(dims(B, order)) do d
        d = DD.maybeshiftlocus(Center(), d)
        # DimPoints are faster if we materialise the dims
        modify(Array, d)
    end
    inpoly = if any(map(d -> DD.order(d) isa Unordered, shifted_dims))
        inpolygon(vec(pts), poly)
    else
        pointbounds = map(shifted_dims) do d
            f = first(d)
            l = last(d)
            f < l ? (f, l) : (l, f)
        end
        # Precalculate vmin, vmax and iyperm permutation vector
        # This is much faster than calling `sortperm` in PolygonInbounds.jl
        vmin = [first.(pointbounds)...]'
        vmax = [last.(pointbounds)...]'
        pts = DimPoints(shifted_dims)
        iyperm = _iyperm(shifted_dims)
        inpolygon(vec(pts), poly; vmin, vmax, iyperm)
    end
    return _inner_fill_polygon!(B, poly, inpoly; order, fill, boundary)
end

function _iyperm(dims::Tuple{<:Dimension,<:Dimension})
    of, ol, os = LA.ordered_firstindex, LA.ordered_lastindex, _order_step
    l1, l2 = map(parent, dims)
    [LinearIndices(size(dims))[i, j] for j in of(l2):os(l2):ol(l2) for i in of(l1):os(l1):ol(l1)]
end
function _iyperm(dims::Tuple{<:Dimension,<:Dimension,<:Dimension})
    of, ol, os = LA.ordered_firstindex, LA.ordered_lastindex, _order_step
    l1, l2, l3 = map(parent, dims)
    [LinearIndices(size(dims))[i, j, k] for k in of(l3):os(l3):ol(l3) for j in of(l2):os(l2):ol(l2) for i in of(l1):os(l1):ol(l1)]
end

_order_step(x) = _order_step(order(x))
_order_step(::ReverseOrdered) = -1
_order_step(::ForwardOrdered) = 1

# split to make a type stability function barrier
function _inner_fill_polygon!(B::AbstractRaster{Bool}, poly, inpoly; order, fill=true, boundary=:center, kw...)
    # Get the array as points
    # Use the first column of the output - the points in the polygon,
    # and reshape to match `A`
    inpolydims = dims(B, order)
    reshaped = Raster(reshape(inpoly, size(inpolydims)), inpolydims)
    for D in DimIndices(B)
        @inbounds if reshaped[D...]
            @inbounds B[D...] = fill
        end
    end
    if boundary === :touches
        # Add the line pixels
        _fill_linestring!(B, poly; order, fill)
    elseif boundary === :inside
        # Remove the line pixels
        _fill_linestring!(B, poly; order, fill=!fill)
    elseif boundary !== :center
        _boundary_error(boundary)
    end
    return B
end

# _fill_point!
# Fill a raster with `fill` where points are inside raster pixels
function _fill_point!(B::AbstractRaster{Bool}, points; order, fill=true, atol=nothing, kw...)
    # Just find which pixels contian the points, and set them to true
    _without_mapped_crs(B) do B1
        for point in points
            selectors = map(dims(B, order), ntuple(i -> i, length(order))) do d, i
                _at_or_contains(d, point[i], atol)
            end
            if hasselection(B, selectors)
                B[selectors...] = fill
            end
        end
    end
    return B
end

# _fill_linestring!
# Fill a raster with `fill` where pixels touch lines in a linestring
function _fill_linestring!(B::AbstractRaster, linestring; order, fill=true, kw...)
    linestring = collect(linestring)
    # Flip the order with a view to keep our alg simple
    forward_ordered_B = reduce(dims(B); init=B) do A, d
        if DD.order(d) isa ReverseOrdered
            A = view(A, rebuild(d, lastindex(d):-1:firstindex(d)))
            set(A, d => reverse(d))
        else
            A
        end
    end
    # For each line segment, burn the line into B with `fill` values
    _fill_line_segments!(forward_ordered_B, linestring, fill, order)
    return B
end

# Separated for a type stability function barrier
function _fill_line_segments!(B::AbstractArray{<:Any,N}, linestring, fill, order) where N
    isfirst = true
    local firstpoint, lastpoint
    for point in linestring
        if isfirst
            isfirst = false
            firstpoint = point
            lastpoint = point
            continue
        end
        if point == firstpoint
            isfirst = true
        end
        line = (
            start=(x=lastpoint[1], y=lastpoint[2]),
            stop=(x=point[1], y=point[2])
        )
        _fill_line!(B, line, fill, order)
        lastpoint = point
    end
end

# fill_line!
# Fill a raster with `fill` where pixels touch a line
# TODO: generalise to 3d and Irregular spacing?
function _fill_line!(A::AbstractRaster, line, fill, order::Tuple{<:Any,<:Any})
    regular = map(dims(A)) do d
        lookup(d) isa Sampled && sampling(d) isa Regular
    end
    all(regular) && throw(ArgumentError( """
            Can only fill lines where dimensions are regular.
            Consider reprojecting the crs, or make an issue in Rasters.jl on github if you need this to work.
            """))

    xd, yd = order
    x_scale = (bounds(A, xd)[2] - bounds(A, xd)[1]) / size(A, xd)
    y_scale = (bounds(A, yd)[2] - bounds(A, yd)[1]) / size(A, yd)
    raw_x_offset = bounds(A, xd)[1]
    raw_y_offset = bounds(A, yd)[1]
    raw_start, raw_stop = line.start, line.stop # Float
    start = (; x=(raw_start.x - raw_x_offset)/x_scale, y=(raw_start.y - raw_y_offset)/y_scale)
    stop = (; x=(raw_stop.x - raw_x_offset)/x_scale, y=(raw_stop.y - raw_y_offset)/y_scale)
    x, y = floor(Int, start.x) + 1, floor(Int, start.y) + 1 # Int

    diff_x = stop.x - start.x
    diff_y = stop.y - start.y
    step_x = signbit(diff_x) * -2 + 1
    step_y = signbit(diff_y) * -2 + 1
    # Ray/Slope calculations
    # Straight distance to the first vertical/horizontal grid boundaries
    xoffset = stop.x > start.x ? (ceil(start.x) - start.x) : (start.x - floor(start.x))
    yoffset = stop.y > start.y ? (ceil(start.y) - start.y) : (start.y - floor(start.y))
    # Angle of ray/slope.
    angle = atan(-diff_y, diff_x)
    # How far to move along the ray to cross the first cell boundary.
    max_x = xoffset / cos(angle)
    max_y = yoffset / sin(angle)
    # How far to move along the ray to move 1 grid cell.
    delta_x = 1.0 / cos(angle)
    delta_y = 1.0 / sin(angle)
    # Travel one grid cell at a time.
    manhattan_distance = floor(Int, abs(floor(start.x) - floor(stop.x)) + abs(floor(start.y) - floor(stop.y)))
    dimconstructors = map(DD.basetypeof, dims(A, order))
    for t in 0:manhattan_distance
        D = map((d, o) -> d(o), dimconstructors, (x, y))
        if checkbounds(Bool, A, D...)
            @inbounds A[D...] = fill
        end
        # Only move in either X or Y coordinates, not both.
        if abs(max_x) <= abs(max_y)
            max_x += delta_x
            x += step_x
        else
            max_y += delta_y
            y += step_y
        end
    end
    return A
end
function _fill_line!(A::AbstractRaster, line, fill, order::Tuple)
    throw(ArgumentError(""""
        Converting a `:line` geometry to raster is currently only implemented for 2d lines.
        Make a Rasters.jl github issue if you need this for more dimensions.
        """
    ))
end

# PolygonInbounds.jl setup

# _flat_nodes
# Convert a geometry/nested vectors to a flat iterator of point nodes for PolygonInbounds
_flat_nodes(A::GI.AbstractGeometry) = _flat_nodes(GI.coordinates(A))
_flat_nodes(A::AbstractVector{<:Union{Missing,<:GI.AbstractGeometry}}) = Iterators.flatten(map(_flat_nodes, A))
function _flat_nodes(A::AbstractVector{<:AbstractVector{<:AbstractVector}})
    Iterators.flatten(map(_flat_nodes, A))
end
_flat_nodes(A::AbstractVector{<:AbstractVector{<:AbstractFloat}}) = A
_flat_nodes(A::AbstractVector{<:Tuple}) = A
_flat_nodes(iter::Base.Iterators.Flatten) = iter


# _to_edges
# Convert a polygon to the `edges` needed by PolygonInbounds
function to_edges(poly)
    edgenum = 0
    edges = Vector{Tuple{Int,Int}}(undef, 0)
    _to_edges!(edges, edgenum, poly)
    return PermutedDimsArray(reinterpret(reshape, Int, edges), (2, 1))
end

# _to_edges!(edges, edgenum, poly)
# fill edges vector with edges from polygon, numbered starting at `edgenum`
function _to_edges!(edges, edgenum, poly::AbstractVector{<:GI.AbstractGeometry})
    foldl(poly; init=edgenum) do n, p
        _to_edges!(edges, n, GI.coordinates(p))
    end
end
function _to_edges!(edges, edgenum, poly::AbstractVector{<:AbstractVector})
    foldl(poly; init=edgenum) do n, p
        _to_edges!(edges, n, p)
    end
end
# Analyse a single polygon
function _to_edges!(
    edges, edgenum, poly::AbstractVector{<:Union{<:NTuple{<:Any,T},<:AbstractVector{T}}}
) where T<:Real
    # Any polygon may actually be made up of multiple sub-polygons,
    # indicated by returning to the first point in the sub-polygon.
    # So we keep track of the first point, and insert closing edges
    # wherever we find it repeated.
    fresh_start = true
    startpoint = first(poly)
    start_edgenum = edgenum
    added_edges = 0
    for (e, point) in enumerate(poly)
        if fresh_start
            # The first edge in the sub-polygon
            start_edgenum = edgenum + e
            startedge = (start_edgenum, start_edgenum + 1)
            push!(edges, startedge)
            startpoint = point
            fresh_start = false
        elseif point == startpoint
            # The closing edge of a sub-polygon 
            closingedge = (edgenum + e, start_edgenum)
            push!(edges, closingedge)
            fresh_start = true
        else
            # A regular edge somewhere in a sub-polygon
            edge = (edgenum + e, edgenum + e + 1)
            push!(edges, edge)
        end
        # Track the total number of edges we have added
        added_edges += 1
    end
    if last(edges)[2] != start_edgenum
        push!(edges, (last(edges)[2], first_edgenum))
        edgecount += 1
    end
    return edgenum + added_edges
end

_flatlength(x::Base.Iterators.Flatten) = x.it
_flatlength(x::AbstractVector) = mapreduce(_flatlength, +, x)
_flatlength(x::AbstractVector{<:Pt}) = length(x)

# _geom_bounds
# Get the bounds of a geometry
function _geom_bounds(geom, order)
    nodes = _flat_nodes(geom)
    bounds = map(ntuple(identity, length(order))) do i
        extrema(p[i] for p in nodes)
    end
    return bounds
end

function _wrapped_geom_bounds(x, geom, order)
    wrapped_bounds = map(rebuild, dims(x, order), _geom_bounds(geom, order))
    return dims(wrapped_bounds, dims(x))
end

# bounds_overlap
# Check if two bounds tuples overlap
function bounds_overlap(a, b)
    axes_cross = map(a, b) do (a_min, a_max), (b_min, b_max)
        if a_max >= b_max
            a_min <= b_max
        else
            a_max >= b_min
        end
    end
    return all(axes_cross)
end

# bbox_might_overlap
# Check if there is a bbox fot the geometry
function bbox_might_overlap(x, geom, order)
    x_bnds = bounds(x, order)
    p_bbox = GI.bbox(geom)
    # If there is no bbox available just act as if it
    # overlaps and let the checks later on sort it out
    p_bbox isa Nothing && return true
    bbox_dims = length(p_bbox) ÷ 2
    p_bnds = [(p_bbox[i], p_bbox[i+bbox_dims]) for i in 1:bbox_dims]

    return bounds_overlap(p_bnds, x_bnds)
end

# _geom_shape
# Get the shape category for a geometry
function _geom_shape(geom)
    typ = GI.geotype(geom)
    if typ in (:Point, :MultiPoint)
        return :point
    elseif typ in (:LineString, :MultiLineString)
        return :line
    elseif typ in (:Polygon, :MultiPolygon)
        return :polygon
    else
        throw(ArgumentError("$typ not known"))
    end
    # TODO: What to do with :GeometryCollection
end

# unwrap_point
# Unwrap a GeoInterface point, to a tuple
function unwrap_point(q::GI.AbstractPoint)
    (q.x, q.y)
end
unwrap_point(q) = q

# Copied from PolygonInbounds, to add extra keyword arguments
function inpoly2(vert, node, edge=zeros(Int);
    atol::T=0.0, rtol::T=NaN, iyperm=nothing, vmin=nothing, vmax=nothing
) where T<:AbstractFloat
    rtol = !isnan(rtol) ? rtol : iszero(atol) ? eps(T)^0.85 : zero(T)
    poly = PolygonInbounds.PolygonMesh(node, edge)
    points = PolygonInbounds.PointsInbound(vert)
    npoints = length(points)

    vmin = isnothing(vmin) ? minimum(points) : vmin
    vmax = isnothing(vmax) ? maximum(points) : vmax
    pmin = minimum(poly)
    pmax = maximum(poly)

    lbar = sum(pmax - pmin)
    tol = max(abs(rtol * lbar), abs(atol))

    ac = PolygonInbounds.areacount(poly)
    stat = ac > 1 ? falses(npoints, 2, ac) : falses(npoints, 2)
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

_shape_error(shape) = throw(ArgumentError("`shape` must be :point, :line or :polygon, got $shape"))
_boundary_error(boundary) = throw(ArgumentError("`boundary` can be :touches, :inside, or :center, got $boundary"))
