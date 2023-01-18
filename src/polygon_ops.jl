
const Poly = AbstractVector{<:Union{NTuple{<:Any,<:Real},AbstractVector{<:Real}}}

const DEFAULT_POINT_ORDER = (XDim, YDim, ZDim)
const DEFAULT_TABLE_DIM_KEYS = (:X, :Y, :Z)

# _fill_geometry!
# Fill a raster with `fill` where it interacts with a geometry.
# This is used in `boolmask` TODO move to mask.jl ?
function fill_geometry!(B::AbstractRaster, data; kw...)
    if Tables.istable(data)
        geomcolname = first(GI.geometrycolumns(data))
        for row in Tables.rows(data)
            geom = Tables.getcolumn(row, geomcolname)
            _fill_geometry!(B, GI.trait(geom), geom; kw...)
        end
    else
        _fill_geometry!(B, GI.trait(data), data; kw...)
    end
end

# This feature filling is simplistic in that it does not use any feature properties.
# This is suitable for masking. See `rasterize` for a version using properties.
_fill_geometry!(B, obj; kw...) = _fill_geometry!(B, GI.trait(obj), obj; kw...)
function _fill_geometry!(B::AbstractRaster, ::GI.AbstractFeatureTrait, feature; kw...)
    _fill_geometry!(B, GI.geometry(feature); kw...)
end
function _fill_geometry!(B::AbstractRaster, ::GI.AbstractFeatureCollectionTrait, fc; kw...)
    _fill_geometry!(B, (GI.geometry(f) for f in GI.getfeature(fc)); kw...)
end
function _fill_geometry!(B::AbstractRaster, ::GI.AbstractGeometryTrait, geom; shape=nothing, verbose=true, kw...)
    shape = shape isa Symbol ? shape : _geom_shape(geom)
    if shape === :point
        _fill_point!(B, geom; shape, kw...)
    elseif shape === :line
        _fill_lines!(B, geom; shape, kw...)
    elseif shape === :polygon
        geomextent = _extent(geom)
        arrayextent = Extents.extent(B, DEFAULT_POINT_ORDER)
        # Only fill if the gemoetry bounding box overlaps the array bounding box
        if !Extents.intersects(geomextent, arrayextent) 
            if verbose
                @info "A geometry was ignored at $geomextent as it was outside of the supplied extent $arrayextent"
            end
            return B
        end
        _fill_polygon!(B, geom; shape, geomextent, kw...)
    else
        throw(ArgumentError("`shape` is $shape, must be `:point`, `:line`, `:polygon` or `nothing`"))
    end
    return B
end
# Treat geoms as an iterator
function _fill_geometry!(B::AbstractRaster, trait::Nothing, geoms; combine=true, kw...)
    if combine
        for geom in geoms
            _fill_geometry!(B, geom; kw...)
        end
    else
        for (i, geom) in enumerate(geoms)
            B1 = view(B, Dim{:geometry}(i))
            _fill_geometry!(B1, geom; kw...)
        end
    end
end

# _fill_polygon!
# Fill a raster with `fill` where pixels are inside a polygon
# `boundary` determines how edges are handled
function _fill_polygon!(B::AbstractRaster, geom;
    fill=true, boundary=:center, geomextent , kw...
)
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
    pts = vec((DimPoints(shifted_dims)))
    inpoly = if any(map(d -> DD.order(d) isa Unordered, shifted_dims))
        inpolygon(pts, geom)
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
        iyperm = _iyperm(shifted_dims)
        inpolygon(pts, geom; vmin, vmax, pmin, pmax, iyperm)
    end
    inpolydims = dims(B, DEFAULT_POINT_ORDER)
    reshaped = Raster(reshape(inpoly, size(inpolydims)), inpolydims)
    return _inner_fill_polygon!(B, geom, reshaped, shifted_dims; fill, boundary)
end

# split to make a type stability function barrier
function _inner_fill_polygon!(B::AbstractRaster, geom, reshaped, shifted_dims; fill::Bool=true, boundary=:center, verbose=true, kw...)
    # Get the array as points
    # Use the first column of the output - the points in the polygon,
    # and reshape to match `A`
    
    # Used to warn if no fill is written for this polygon
    any_in_polygon = false
    n_on_line = false

    # TODO: This takes a while, it could be faster to use
    # modified PolygonInbounds output directly rather than 
    # a reshaped view of the first column?
    for D in DimIndices(B)
        @inbounds is_in_polygon = reshaped[D...]
        if is_in_polygon 
            any_in_polygon = true
            @inbounds B[D...] = fill
        end
    end
    if boundary === :touches
        B1 = setdims(B, shifted_dims)
        # Add the line pixels
        n_on_line = _fill_lines!(B1, geom; fill)
    elseif boundary === :inside
        B1 = setdims(B, shifted_dims)
        # Remove the line pixels
        n_on_line = _fill_lines!(B1, geom; fill=!fill)
    elseif boundary !== :center
        throw(ArgumentError("`boundary` can be :touches, :inside, or :center, got $boundary"))
    end
    if verbose
        (any_in_polygon | n_on_line > 0) || @warn "$n polygons were not filled as they did not cross a pixel center. Consider using `boundary=:touches`, or `verbose=false` to hide these warning"
    end
    return B
end

struct IYPerm{D<:Tuple,R} <: AbstractVector{Int}
    dims::D
    ranges::R
end
function IYPerm(dims::D) where D
    ranges = map(dims) do d
        l = parent(d)
        LA.ordered_firstindex(l):_order_step(l):LA.ordered_lastindex(l)
    end
    IYPerm{D,typeof(ranges)}(dims, ranges)
end
DD.dims(yp::IYPerm) = yp.dims

_order_step(x) = _order_step(order(x))
_order_step(::ReverseOrdered) = -1
_order_step(::ForwardOrdered) = 1

Base.@propagate_inbounds function Base.getindex(yp::IYPerm, i::Int)
    ci = Tuple(CartesianIndices(size(dims(yp)))[i])
    I = map(getindex, yp.ranges, ci)
    return LinearIndices(size(dims(yp)))[I...]
end
Base.size(yp::IYPerm) = (prod(map(length, dims(yp))),)

_iyperm(dims::Tuple) = IYPerm(dims)

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

# _fill_lines!
# Fill a raster with `fill` where pixels touch lines in a geom
# Separated for a type stability function barrier
function _fill_lines!(B::AbstractRaster, geom; fill=true, kw...)
    # Make sure dims have `Center` locus
    centered_dims = map(dims(B, DEFAULT_POINT_ORDER)) do d
        d = DD.maybeshiftlocus(Center(), d)
    end
    centered_B = setdims(B, centered_dims)
    # Flip the order with a view to keep our alg simple
    forward_ordered_B = reduce(dims(centered_B); init=centered_B) do A, d
        if DD.order(d) isa ReverseOrdered
            A = view(A, rebuild(d, lastindex(d):-1:firstindex(d)))
            set(A, d => reverse(d))
        else
            A
        end
    end
    return _fill_lines!(forward_ordered_B, geom, fill)
end

_fill_lines!(forward_ordered_B, geom, fill) =
    _fill_lines!(forward_ordered_B, GI.geomtrait(geom), geom, fill)
function _fill_lines!(B::AbstractArray, ::Union{GI.MultiLineStringTrait}, geom, fill)
    n_on_line = 0
    for linestring in GI.getlinestring(geom)
        n_on_line += _fill_lines!(B, linestring, fill)
    end
    return n_on_line
end
function _fill_lines!(
    B::AbstractArray, ::Union{GI.MultiPolygonTrait,GI.PolygonTrait}, geom, fill
)
    n_on_line = 0
    for ring in GI.getring(geom)
        n_on_line += _fill_lines!(B, ring, fill)
    end
    return n_on_line
end
function _fill_lines!(
    B::AbstractArray, ::Union{GI.LineStringTrait,GI.LinearRingTrait}, linestring, fill
)
    isfirst = true
    local firstpoint, lastpoint
    n_on_line = 0
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
        n_on_line += _burn_line!(B, line, fill)
        lastpoint = point
    end
    return n_on_line
end

# _burn_line!
#
# Line-burning algorithm
# Burns a line into a raster with `fill` value where pixels touch a line
#
# TODO: generalise to Irregular spans?
function _burn_line!(A::AbstractRaster, line, fill)
    xd, yd = dims(A, DEFAULT_POINT_ORDER)
    regular = map((xd, yd)) do d
        lookup(d) isa AbstractSampled && span(d) isa Regular
    end
    msg = """
        Can only fill lines where dimensions have `Regular` lookups.
        Consider using `boundary=:center`, reprojecting the crs,
        or make an issue in Rasters.jl on github if you need this to work.
        """
    all(regular) || throw(ArgumentError(msg))

    x_scale = abs(step(span(A, X)))
    y_scale = abs(step(span(A, Y)))
    raw_x_offset = xd[1] - x_scale / 2
    raw_y_offset = yd[1] - x_scale / 2
    # Converted lookup to array axis values
    start = (x=(line.start.x - raw_x_offset)/x_scale, y=(line.start.y - raw_y_offset)/y_scale)
    stop = (x=(line.stop.x - raw_x_offset)/x_scale, y=(line.stop.y - raw_y_offset)/y_scale)
    x, y = ceil(Int, start.x), ceil(Int, start.y) # Int
    @show x, y

    diff_x = stop.x - start.x
    diff_y = stop.y - start.y
    step_x = signbit(diff_x) * -2 + 1
    step_y = signbit(diff_y) * -2 + 1
    # Ray/Slope calculations
    # Straight distance to the first vertical/horizontal grid boundaries
    xoffset = stop.x > start.x ? (ceil(start.x) - start.x) : (start.x - floor(start.x))
    yoffset = stop.y > start.y ? (ceil(start.y) - start.y) : (start.y - floor(start.y))
    @show xoffset yoffset
    # Angle of ray/slope.
    # max: How far to move along the ray to cross the first cell boundary.
    # delta: How far to move along the ray to move 1 grid cell.
    # Our precision requirements are very low given pixels per
    # axis is usually some number of thousands
    ang = @fastmath atan(-diff_y, diff_x)
    cs = @fastmath cos(ang)
    si = @fastmath sin(ang)
    delta_x, max_x = if isapprox(cs, zero(cs); atol=1e-10)
        @show Inf
        -Inf, Inf
    else
        1.0 / cs, xoffset / cs
    end
    delta_y, max_y = if isapprox(si, zero(si); atol=1e-10)
        @show Inf
        -Inf, Inf
    else
        1.0 / si, yoffset / si
    end
    # For arbitrary dimension indexing
    dimconstructors = map(DD.basetypeof, (xd, yd))
    # Count how many exactly hit lines
    n_on_line = 0
    # Travel one grid cell at a time.
    manhattan_distance = floor(Int, abs(ceil(start.x) - floor(stop.x)) + abs(ceil(start.y) - floor(stop.y)))
    @show manhattan_distance
    for _ in 0:manhattan_distance
        D = map((d, o) -> d(o), dimconstructors, (x, y))
        if checkbounds(Bool, A, D...)
            n_on_line += 1
            val = fill isa Function ? fill(A[D...]) : fill
            @inbounds A[D...] = val
        end
        # Only move in either X or Y coordinates, not both.
        if abs(max_x) < abs(max_y)
            max_x += delta_x
            x += step_x
        else
            max_y += delta_y
            y += step_y
        end
    end
    return n_on_line
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
    for (n, point) in enumerate(GI.getpoint(geom))
        i = lastnode + n
        if n == npoints
            # The closing edge of a sub-polygon
            edges[i, 1] = i
            edges[i, 2] = lastnode + 1
        else
            # A regular edge somewhere in a sub-polygon
            edges[i, 1] = i
            edges[i, 2] = i + 1
        end
        nodes[i, 1] = GI.x(point)
        nodes[i, 2] = GI.y(point)
    end
    return lastnode + npoints
end

# _extent
# Get the bounds of a geometry
_extent(geom) = _extent(GI.trait(geom), geom)
function _extent(::Nothing, data::AbstractVector)
    reduce(data; init=_extent(first(data))) do ext, geom
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
            # TODO: test this branch
            # Table of points with dimension columns
            data = reduce(DEFAULT_TABLE_DIM_KEYS; init=(;)) do acc, key
                if key in Tables.columnnames(cols)
                    merge(acc, (; key=extrema(cols[key])))
                else
                    acc
                end
            end
            return Extent(data)
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
    # TODO: this is not actually used.
    rasterextent = Extents.extent(x, DEFAULT_POINT_ORDER)
    geomextent = GI.extent(geom)
    if isnothing(geomextent)
        return true
    else
        return Extents.intersects(geomextent, rasterextent)
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
