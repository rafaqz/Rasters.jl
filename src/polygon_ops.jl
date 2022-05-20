const Pt{T<:Real} = Union{AbstractVector{T},NTuple{<:Any,T}}
const Poly = AbstractVector{<:Union{NTuple{<:Any,<:Real},AbstractVector{<:Real}}}

const DEFAULT_POINT_ORDER = (XDim, YDim)

# _fill_geometry!
# Fill a raster with `fill` where it interacts with a geometry.
function _fill_geometry!(B::AbstractRaster, geom; order, kw...)
    if isfeature(geom)
        bbox_might_overlap(B, geom, order) || return B
        _fill_geometry!(B, GI.geometry(geom), geom; order, kw...)
    elseif isgeom(geom)
        bbox_might_overlap(B, geom, order) || return B
        _fill_geometry!(B, geom; order, shape=_geom_shape(geom), kw...)
    end
end
function _fill_geometry!(B::AbstractRaster, geoms::AbstractVector; shape=:polygon, order, kw...)
    if isfeature(first(geoms))
        bbox_might_overlap(B, geom, order) || return B
        for feature in geom
            _fill_geometry!(B, feature; order, kw...)
        end
        return B
    elseif isgeometry(first(geoms))
        any(map(g -> bbox_might_overlap(B, g, order), geoms)) || return B
        if eltype(geoms) == typeof(first(geoms))
            shape = _geom_shape(first(geom))
            for geom in geoms
                _fill_geometry!(B, geom; order, shape, kw...)
            end
        else
            _fill_mixed_geometries!(B, geoms; order, kw...)
        end
        return B
    else
        gbounds = _geom_bounds(geoms, order)
        abounds = bounds(dims(B, order)) # Only mask if the gemoetry bounding box overlaps the array bounding box
        bounds_overlap(gbounds, abounds) || return B
        if shape === :polygon
            _fill_polygon!(B, geoms; polybounds=gbounds, order, kw...)
        elseif shape === :line
            _fill_linestring!(B, geoms; order, kw...)
        elseif shape === :point
            _fill_points!(B, geoms; order, kw...)
        else
            _shape_error(shape)
        end
        return B
    end
end

# Multiple geometry types. Rasterize by groups of types separately
@noinline function _fill_mixed_geometries!(B, geoms; order, kw...)
    geomgroups = Dict{Type,Vector}()
    for geom in geoms
        key = typeof(geom)
        if haskey(geomgroups, key)
            push!(geomgroups[key], geom)  
        else
            geomgroups[key] = [geom]
        end
    end
    for group in values(geomgroups)
        _fill_geometry!(B, GI.coordinates.(group); order, shape=_geom_shape(first(group)), kw...)
    end
end

# _fill_polygon!
# Fill a raster with `fill` where pixels are inside a polygon
# `boundary` determines how edges are handled
function _fill_polygon!(B::AbstractRaster, poly; polybounds, order, fill=true, boundary=:center, kw...)
    # TODO take a view of B for the polygon
    # We need a tuple of all the dims in `order`
    # We also need the index locus to be the center so we are
    # only selecting cells more than half inside the polygon
    shifted_dims = map(dims(B, order)) do d
        d = DD.maybeshiftlocus(Center(), d)
        # DimPoints are faster if we materialise the dims
        modify(Array, d)
    end
    pts = DimPoints(shifted_dims)
    points = getpoint(polygon)
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
        pmin = [first.(polybounds)...]'
        pmax = [last.(polybounds)...]'
        pts = DimPoints(shifted_dims)
        iyperm = _iyperm(shifted_dims)
        inpolygon(vec(pts), poly; vmin, vmax, pmin, pmax, iyperm)
    end
    inpolydims = dims(B, order)
    reshaped = Raster(reshape(inpoly, size(inpolydims)), inpolydims)
    return _inner_fill_polygon!(B, poly, inpoly, reshaped; order, fill, boundary)
end

# split to make a type stability function barrier
function _inner_fill_polygon!(B::AbstractRaster, poly, inpoly, reshaped; order, fill=true, boundary=:center, kw...)
    # Get the array as points
    # Use the first column of the output - the points in the polygon,
    # and reshape to match `A`
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

function _iyperm(dims::Tuple{<:Dimension,<:Dimension})
    a1, a2 = map(dims) do d
        l = parent(d)
        LA.ordered_firstindex(l):_order_step(l):LA.ordered_lastindex(l)
    end
    iyperm = Array{Int}(undef, length(a1) * length(a2))
    lis = (LinearIndices(size(dims))[i, j] for j in a2 for i in a1)
    for (i, li) in enumerate(lis)
        iyperm[i] = li
    end
    return iyperm
end
function _iyperm(dims::Tuple{<:Dimension,<:Dimension,<:Dimension})
    a1, a2, a3 = map(dims) do d
        l = parent(d)
        LA.ordered_firstindex(l):_order_step(l):LA.ordered_lastindex(l)
    end
    iyperm = Array{Int}(undef, length(a1) * length(a2) * length(a3))
    lis = (LinearIndices(size(dims))[i, j, k] for k in a3 for j in a2 for i in a1)
    for (i, li) in enumerate(lis)
        Iyperm[i] = li
    end
    return iyperm
end

_order_step(x) = _order_step(order(x))
_order_step(::ReverseOrdered) = -1
_order_step(::ForwardOrdered) = 1

# _fill_point!
# Fill a raster with `fill` where points are inside raster pixels
function _fill_points!(B::AbstractRaster, geom; kw...)
    # Just find which pixels contian the points, and set them to true
    _without_mapped_crs(B) do B1
        for point in getpoint(geom)
            _fill_point!(B, point; kw...)
        end
    end
    return B
end

function _fill_point!(B::AbstractRaster, point; 
    order, fill=true, atol=nothing, kw...
)
    selectors = map(dims(B, order), ntuple(i -> i, length(order))) do d, i
        _at_or_contains(d, getcoord(point, i], atol)
    end
    if hasselection(B, selectors)
        B[selectors...] = fill
    end
end

# _fill_linestring!
# Fill a raster with `fill` where pixels touch lines in a linestring
# Separated for a type stability function barrier
function _fill_linestring!(B::AbstractRaster, linestring::AbstractVector{<:AbstractVector{<:AbstractVector}}, args...; kw...)
    for ls in linestring
        _fill_linestring!(B, ls, args...; kw...)
    end
end
function _fill_linestring!(B::AbstractRaster, linestring::Poly; order, fill=true, kw...)
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

function _fill_line_segments!(B::AbstractArray, linestring::Poly, fill, order::Tuple{<:Dimension,<:Dimension})
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
function _fill_line!(A::AbstractRaster, line, fill, order::Tuple{<:Dimension,<:Dimension})
    regular = map(dims(A, order)) do d
        lookup(d) isa AbstractSampled && span(d) isa Regular
    end
    msg = """
        Can only fill lines where dimensions are regular.
        Consider using `boundary=center` reprojecting the crs,
        or make an issue in Rasters.jl on github if you need this to work.
        """
    all(regular) || throw(ArgumentError(msg))

    xd, yd = order
    x_scale = abs(step(span(A, X)))
    y_scale = abs(step(span(A, Y)))
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
    # max: How far to move along the ray to cross the first cell boundary.
    # delta: How far to move along the ray to move 1 grid cell.
    cs = cos(angle)
    si = sin(angle)
    max_x, delta_x = if isapprox(cs, zero(cs); atol=1e-10) 
        -Inf, Inf
    else
        1.0 / cs, xoffset / cs
    end
    max_y, delta_y = if isapprox(si, zero(si); atol=1e-10)
        -Inf, Inf
    else
        1.0 / si, yoffset / si
    end
    # Travel one grid cell at a time.
    manhattan_distance = floor(Int, abs(floor(start.x) - floor(stop.x)) + abs(floor(start.y) - floor(stop.y)))
    # For arbitrary dimension indexing
    dimconstructors = map(DD.basetypeof, dims(A, order))
    for t in 0:manhattan_distance
        D = map((d, o) -> d(o), dimconstructors, (x, y))
        if checkbounds(Bool, A, D...)
            if fill isa Function 
                @inbounds A[D...] = fill(A[D...])
            else
                @inbounds A[D...] = fill
            end
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
function _fill_line!(A::AbstractRaster, line, fill, order::Tuple{Vararg{<:Dimension}})
    msg = """"
        Converting a `:line` geometry to raster is currently only implemented for 2d lines.
        Make a Rasters.jl github issue if you need this for more dimensions.
        """
    throw(ArgumentError(msg))
end

# PolygonInbounds.jl setup

# _flat_nodes
# Convert a geometry/nested vectors to a flat iterator of point nodes for PolygonInbounds
function _flat_nodes(A::AbstractVector{T})
    if GI.geomtrait(first(A)) isa AbstractPointTrait
        A
    else
        Iterators.flatten(map(_flat_nodes, A))
    end
end
_flat_nodes(iter::Base.Iterators.Flatten) = iter

# _to_edges
# Convert a polygon to the `edges` needed by PolygonInbounds
to_edges_and_nodes(poly) = to_edges_and_nodes(GI.geomtrait(poly), poly)
function to_edges_and_nodes(::Nothing, poly)
    p1 = first(_flat_nodes(poly))
    edges = Vector{Tuple{Int,Int}}(undef, 0)
    nodes = Vector{typeof(p1)}(undef, 0)
    nodenum = 0
    _to_edges!(edges, nodes, nodenum, poly)
    edges = PermutedDimsArray(reinterpret(reshape, Int, edges), (2, 1))
    lastx = 0
    for (i, x) in enumerate(view(edges, :, 1))
        x == lastx + 1 || error("edge dup at $x $lastx $i")
        lastx = x
    end
    return edges, nodes
end
function to_edges_and_nodes(::Union{AbstractPolygonTrait,AbstractMultiPolygonTrait}, poly)
    n = npoint(poly)
    edges = Vector{Tuple{Int,Int}}(undef, n)
    nodes = Vector{typeof(p1)}(undef, n)
    nodenum = 0
    for ring in getring(poly)
        nodenum = _to_edges!(edges, nodes, nodenum, poly)
    end
end

# _to_edges!(edges, edgenum, poly)
# fill edges vector with edges from polygon, numbered starting at `edgenum`
function _to_edges!(edges, nodes, pointnum, poly::AbstractVector)
    foldl(poly; init=pointnum) do n, p
        _to_edges!(edges, nodes, n, p)
    end
end
# Analyse a single polygon
function _to_edges!(edges, nodes, pointnum, poly)
    fresh_start = true
    startpoint = getpoint(poly, 1)
    start_pointnum = pointnum
    added_nodes = 0
    for (n, point) in enumerate(poly)
        if fresh_start
            # The first edge in the sub-polygon
            start_pointnum = pointnum + n
            startedge = (start_pointnum, start_pointnum + 1)
            push!(edges, startedge)
            startpoint = point
            fresh_start = false
        elseif point == startpoint || n == length(poly)
            # The closing edge of a sub-polygon 
            closingedge = (pointnum + n, start_pointnum)
            push!(edges, closingedge)
            fresh_start = true
        else
            # A regular edge somewhere in a sub-polygon
            edge = (pointnum + n, pointnum + n + 1)
            push!(edges, edge)
        end
        # Track the total number of nodes we have added
        added_nodes = n
    end
    append!(nodes, poly)
    nextpoint = pointnum + added_nodes
    return nextpoint
end

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
# Check if there is a bbox for the geometry
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
_geom_shape(geom) = _geom_shape(GI.geomtrait(geom), geom)
_geom_shape(geom::Union{PointTrait,MultiPointTrait}) = :point
_geom_shape(geom::Union{LineStringTrait,MultiLineStringTrait}) = :line 
_geom_shape(geom::Union{LinearRingTrait,PolygonTrait,MultiPolygonTrait}) = :polygon
_geom_shape(trait, geom) = throw(ArgumentError("Geometry trait $trait not handled by Rasters.jl"))

# Copied from PolygonInbounds, to add extra keyword arguments
# PR to include these when this has solidied
function inpoly2(vert, node, edge=zeros(Int);
    atol::T=0.0, rtol::T=NaN, iyperm=nothing, 
    vmin=nothing, vmax=nothing, pmin=nothing, pmax=nothing
) where T<:AbstractFloat
    rtol = !isnan(rtol) ? rtol : iszero(atol) ? eps(T)^0.85 : zero(T)
    poly = PolygonInbounds.PolygonMesh(node, edge)
    points = PolygonInbounds.PointsInbound(vert)
    npoints = length(points)

    vmin = isnothing(vmin) ? minimum(points) : vmin
    vmax = isnothing(vmax) ? maximum(points) : vmax
    pmin = isnothing(pmin) ? minimum(poly) : pmin
    pmax = isnothing(pmax) ? maximum(poly) : pmax

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
