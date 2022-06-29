
const Poly = AbstractVector{<:Union{NTuple{<:Any,<:Real},AbstractVector{<:Real}}}

const DEFAULT_POINT_ORDER = (XDim, YDim, ZDim)
const DEFAULT_TABLE_DIM_KEYS = (:X, :Y, :Z)

# _fill_geometry!
# Fill a raster with `fill` where it interacts with a geometry.
function _fill_geometry!(B::AbstractRaster, geom; kw...)
    _fill_geometry!(B, GI.trait(geom), geom; kw...)
end
function _fill_geometry!(B::AbstractRaster, ::GI.AbstractFeatureTrait, geom; kw...)
    return _fill_geometry!(B, GI.geometry(geom), geom; kw...)
end
function _fill_geometry!(B::AbstractRaster, ::GI.AbstractFeatureCollectionTrait, fc; kw...)
    for feature in GI.getfeature(fc)
        _fill_geometry!(B, GI.geometry(feature); kw...)
    end
end
function _fill_geometry!(B::AbstractRaster, ::GI.AbstractGeometryTrait, geom; shape=nothing, kw...)
    shape = shape isa Symbol ? shape : _geom_shape(geom)
    if shape === :point
        _fill_point!(B, geom; shape, kw...)
    elseif shape === :line
        _fill_linestring!(B, geom; shape, kw...)
    elseif shape === :polygon
        geomextent = _extent(geom)
        arrayextent = Extents.extent(B, DEFAULT_POINT_ORDER) 
        # Only fill if the gemoetry bounding box overlaps the array bounding box
        Extents.intersects(geomextent, arrayextent) || return B
        _fill_polygon!(B, geom; shape, geomextent, kw...)
    else
        throw(ArgumentError("`shape` is $shape, must be `:point`, `:line`, `:polygon` or `nothing`"))
    end
    return B
end
function _fill_geometry!(B::AbstractRaster, trait::Nothing, geoms::AbstractVector; kw...)
    for geom in geoms
        _fill_geometry!(B, geom; kw...)
    end
end
function _fill_geometry!(B::AbstractRaster, trait::Nothing, data; kw...)
    if Tables.istable(data)
        geomcolname = first(GI.geometrycolumns(data))
        for row in Tables.rows(data)
            geom = Tables.getcolumn(row, geomcolname)
            _fill_geometry!(B, geom; kw...)
        end
    else
        throw(ArgumentError("Unknown geometry object $(typeof(data))"))
    end
end

# _fill_polygon!
# Fill a raster with `fill` where pixels are inside a polygon
# `boundary` determines how edges are handled
function _fill_polygon!(B::AbstractRaster, geom; geomextent, fill=true, boundary=:center, kw...)
    # Subset to the area the geom covers
    # TODO take a view of B for the polygon
    # We need a tuple of all the dims in `order`
    # We also need the index locus to be the center so we are
    # only selecting cells more than half inside the polygon
    shifted_dims = map(dims(B, DEFAULT_POINT_ORDER)) do d
        d = DD.maybeshiftlocus(Center(), d)
        # DimPoints are faster if we materialise the dims
        modify(Array, d)
    end
    pts = DimPoints(shifted_dims)
    points = GI.getpoint(geom)
    inpoly = if any(map(d -> DD.order(d) isa Unordered, shifted_dims))
        inpolygon(vec(pts), geom)
    else
        pointbounds = map(shifted_dims) do d
            f = first(d)
            l = last(d)
            f < l ? (f, l) : (l, f)
        end
        # Precalculate vmin, vmax and iyperm permutation vector
        # This is much faster than calling `sortperm` in PolygonInbounds.jl
        vmin = [map(first, pointbounds)...]'
        vmax = [map(last, pointbounds)...]'
        pmin = [map(first, geomextent)...]'
        pmax = [map(last, geomextent)...]'
        pts = DimPoints(shifted_dims)
        iyperm = _iyperm(shifted_dims)
        inpolygon(vec(pts), geom; vmin, vmax, pmin, pmax, iyperm)
    end
    inpolydims = dims(B, DEFAULT_POINT_ORDER)
    reshaped = Raster(reshape(inpoly, size(inpolydims)), inpolydims)
    return _inner_fill_polygon!(B, geom, inpoly, reshaped; fill, boundary)
end

# split to make a type stability function barrier
function _inner_fill_polygon!(B::AbstractRaster, geom, inpoly, reshaped; fill=true, boundary=:center, kw...)
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
        _fill_linestring!(B, geom; fill)
    elseif boundary === :inside
        # Remove the line pixels
        _fill_linestring!(B, geom; fill=!fill)
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
@noinline _fill_point!(x::RasterStackOrArray, geom; kw...) = _fill_point!(x, GI.geomtrait(geom), geom; kw...)
@noinline function _fill_point!(x::RasterStackOrArray, ::GI.AbstractGeometryTrait, geom; kw...)
    # Just find which pixels contain the points, and set them to true
    _without_mapped_crs(x) do x1
        for point in GI.getpoint(geom)
            _fill_point!(x, point; kw...)
        end
    end
    return x
end
@noinline function _fill_point!(x::RasterStackOrArray, ::GI.AbstractPointTrait, point;
    fill=true, atol=nothing, kw...
)
    selectors = map(dims(x, DEFAULT_POINT_ORDER)) do d
        _at_or_contains(d, _dimcoord(d, point), atol)
    end
    # TODO make a check in dimensionaldata that returns the index if it is inbounds
    if hasselection(x, selectors)
        I = dims2indices(x, selectors)
        _fill_index!(x, fill, I)
    end
    return x
end

# Fill Int indices directly
function _fill_index!(st::AbstractRasterStack, fill, I)
    foreach(st, fill) do A, f
        _fill_index!(A, f, I)
    end
end
_fill_index!(A::AbstractRaster, fill, I::NTuple{<:Any,Int}) = A[I...] = fill
_fill_index!(A::AbstractRaster, fill::Function, I::NTuple{<:Any,Int}) =
    A[I...] = fill(A[I...])
# Handle filling over arbitrary dimensions.
function _fill_index!(A::AbstractRaster, fill, I)
    v = view(A, I...)
    for i in eachindex(v)
        val = fill isa Function ? fill(v[i]) : fill
        v[i] = val
    end
end

# _fill_linestring!
# Fill a raster with `fill` where pixels touch lines in a linestring
# Separated for a type stability function barrier
function _fill_linestring!(B::AbstractRaster, linestring; fill=true, kw...)
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
    _fill_line_segments!(forward_ordered_B, linestring, fill)
    return B
end

function _fill_line_segments!(B::AbstractArray, linestring, fill)
    isfirst = true
    local firstpoint, lastpoint
    for point in GI.getpoint(linestring)
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
            start=(x=GI.x(lastpoint), y=GI.y(lastpoint)),
            stop=(x=GI.x(point), y=GI.y(point)),
        )
        _fill_line!(B, line, fill)
        lastpoint = point
    end
end

# fill_line!
# Fill a raster with `fill` where pixels touch a line
# TODO: generalise to Irregular spans?
function _fill_line!(A::AbstractRaster, line, fill)
    xd, yd = dims(A, DEFAULT_POINT_ORDER)
    regular = map((xd, yd)) do d
        lookup(d) isa AbstractSampled && span(d) isa Regular
    end
    msg = """
        Can only fill lines where dimensions have `Regular` lookups.
        Consider using `boundary=center`, reprojecting the crs,
        or make an issue in Rasters.jl on github if you need this to work.
        """
    all(regular) || throw(ArgumentError(msg))

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
    angle = @fastmath atan(-diff_y, diff_x)
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
    dimconstructors = map(DD.basetypeof, (xd, yd))
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

# _to_edges
# Convert a polygon to the `edges` needed by PolygonInbounds
to_edges_and_nodes(geom) = to_edges_and_nodes(GI.geomtrait(geom), geom)
function to_edges_and_nodes(
    tr::Union{GI.LinearRingTrait,GI.AbstractPolygonTrait,GI.AbstractMultiPolygonTrait}, geom
)
    n = GI.npoint(geom)
    edges = Matrix{Int}(undef, n, 2)
    nodes = Matrix{Float64}(undef, n, 2)
    nodenum = 0
    if tr isa GI.LinearRingTrait
        to_edges_and_nodes!(edges, nodes, nodenum, geom)
    else
        for ring in GI.getring(geom)
            nodenum = to_edges_and_nodes!(edges, nodes, nodenum, ring)
        end
    end
    return edges, nodes
end
# _to_edges!(edges, edgenum, geom)
# Analyse a single geometry
function to_edges_and_nodes!(edges, nodes, lastnode, geom)
    npoints = GI.npoint(geom)
    for (n, point) in enumerate(GI.getgeom(geom))
        i = lastnode + n
        if n == npoints
            # The closing edge of a sub-polygon 
            closingedge = (i, lastnode + 1)
            edges[i, :] .= closingedge
        else
            # A regular edge somewhere in a sub-polygon
            edge = (i, i + 1)
            edges[i, :] .= edge
        end
        nodes[i, :] .= (GI.x(point), GI.y(point))
    end
    return lastnode + npoints
end

# _extent
# Get the bounds of a geometry
_extent(geom) = _extent(GI.trait(geom), geom)
function _extent(::Nothing, data::AbstractVector)
    reduce(geom; init=_extent(first(data))) do ext, geom
        Extents.union(ext, _extent(geom))
    end
end
_extent(::Nothing, data::RasterStackOrArray) = Extents.extent(data)
function _extent(::Nothing, data)
    if Tables.istable(data)
        geomcolname = first(GI.geometrycolumns(data))
        cols = Tables.columns(data)
        if geomcolname in Tables.columnnames(cols)
            # Table of geometries
            geoms = Tables.getcolumn(data, geomcolname)
            return reduce(geoms; init=_extent(first(geoms))) do ext, geom
                Extents.union(ext, _extent(geom))
            end
        else
            # Table of points with dimension columns
            reduce(DEFAULT_TABLE_DIM_KEYS; init=(;)) do acc, key
                if key in Tables.columnnames(cols) 
                    merge(acc, (; key=extrema(columns)))
                else
                    acc
                end
            end
        end
    else
        return Extents.extent(data)
    end
end
_extent(::GI.AbstractPointTrait, geom) = Extents.Extent(X=GI.x(geom), Y=GI.y(geom))
function _extent(::GI.AbstractTrait, geom)
    geomextent = GI.extent(geom)
    if isnothing(geomextent)
        points = GI.getpoint(geom)
        xbounds = extrema(GI.x(p) for p in points)
        ybounds = extrema(GI.y(p) for p in points)
        return Extents.Extent(X=xbounds, Y=ybounds)
    else
        return geomextent
    end
end
_extent(::GI.AbstractFeatureTrait, feature) = _extent(GI.geometry(feature))
function _extent(::GI.AbstractFeatureCollectionTrait, features) 
    features = GI.getfeature(features)
    reduce(features; init=_extent(first(features))) do acc, f
        Extents.union(acc, _extent(f))
    end
end

# extent_may_intersect
# Check if there is an extent for the geometry
function extent_may_intersect(x, geom)
    rasterextent = Extents.extent(x, DEFAULT_POINT_ORDER)
    geomextent = GI.extent(geom)
    if ext isa Nothing 
        return true
    else
        return Extents.intersects(geomextent, rasterext)
    end
end

# _dimcoord
# Get the GeoInterface coord from a point for a specific Dimension
_dimcoord(::XDim, point) = GI.x(point)
_dimcoord(::YDim, point) = GI.y(point)
_dimcoord(::ZDim, point) = GI.z(point)

# _geom_shape
# Get the shape category for a geometry
@inline _geom_shape(geom) = _geom_shape(GI.geomtrait(geom))
@inline _geom_shape(geom::Union{<:GI.PointTrait,<:GI.MultiPointTrait}) = :point
@inline _geom_shape(geom::Union{<:GI.LineStringTrait,<:GI.MultiLineStringTrait}) = :line 
@inline _geom_shape(geom::Union{<:GI.LinearRingTrait,<:GI.PolygonTrait,<:GI.MultiPolygonTrait}) = :polygon
# _geom_shape(trait, geom) = throw(ArgumentError("Geometry trait $trait not handled by Rasters.jl"))

_shape_error(shape) = throw(ArgumentError("`shape` must be :point, :line or :polygon, got $shape"))
_boundary_error(boundary) = throw(ArgumentError("`boundary` can be :touches, :inside, or :center, got $boundary"))
