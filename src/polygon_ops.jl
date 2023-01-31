
const Poly = AbstractVector{<:Union{NTuple{<:Any,<:Real},AbstractVector{<:Real}}}

const DEFAULT_POINT_ORDER = (X(), Y())
const DEFAULT_TABLE_DIM_KEYS = (:X, :Y)

# _burn_geometry!
# Fill a raster with `fill` where it interacts with a geometry.
# This is used in `boolmask` TODO move to mask.jl ?
# 
# _istable keyword is a hack so we know not to pay the
# price of calling `istable` which calls `hasmethod`
function burn_geometry!(B::AbstractRaster, data; kw...)
    if Tables.istable(typeof(data))
        geomcolname = first(GI.geometrycolumns(data))
        for row in Tables.rows(data)
            geom = Tables.getcolumn(row, geomcolname)
            _burn_geometry!(B, GI.trait(geom), geom; kw...)
        end
    else
        _burn_geometry!(B, GI.trait(data), data; kw...)
    end
    return B
end

# This feature filling is simplistic in that it does not use any feature properties.
# This is suitable for masking. See `rasterize` for a version using properties.
_burn_geometry!(B, obj; kw...) = _burn_geometry!(B, GI.trait(obj), obj; kw...)
function _burn_geometry!(B::AbstractRaster, ::GI.AbstractFeatureTrait, feature; kw...)
    _burn_geometry!(B, GI.geometry(feature); kw...)
end
function _burn_geometry!(B::AbstractRaster, ::GI.AbstractFeatureCollectionTrait, fc; kw...)
    _burn_geometry!(B, (GI.geometry(f) for f in GI.getfeature(fc)); kw...)
end
function _burn_geometry!(B::AbstractRaster, ::GI.AbstractGeometryTrait, geom; shape=nothing, _buffer=nothing, verbose=false, kw...)
    # Use the specified shape or detect it
    shape = shape isa Symbol ? shape : _geom_shape(geom)
    if shape === :point
        _fill_point!(B, geom; fill=true, shape, kw...)
    elseif shape === :line
        _burn_lines!(B, geom; shape, kw...)
    elseif shape === :polygon
        # Get the extents of the geometry and array
        geomextent = _extent(geom)
        arrayextent = Extents.extent(B, DEFAULT_POINT_ORDER)
        # Only fill if the gemoetry bounding box overlaps the array bounding box
        if !Extents.intersects(geomextent, arrayextent) 
            if verbose
                @info "A geometry was ignored at $geomextent as it was outside of the supplied extent $arrayextent"
            end
            return B
        end
        # Take a view of the geometry extent
        B1 = view(B, Touches(geomextent))
        if isnothing(_buffer)
            # Burn straight into the main array (for one object rasterization)
            _burn_polygon!(B1, geom; shape, geomextent, kw...)
        else
            # Also take a view of the buffer
            buf1 = view(_buffer, Touches(geomextent))
            # Reset the buffer for this area
            buf1 .= false
            # Burn the polygon into the buffer
            _burn_polygon!(buf1, geom; shape, geomextent, kw...)
            # btwise or any burned pixels into the main array
            B1 .|= buf1
        end
    else
        throw(ArgumentError("`shape` is $shape, must be `:point`, `:line`, `:polygon` or `nothing`"))
    end
    return B
end
# Treat geoms as an iterator
function _burn_geometry!(B::AbstractRaster, trait::Nothing, geoms; combine=true, kw...)
    if combine
        # Define a buffer to write to, we cant burn lines and polygons
        # into the same buffer - they reverse each other.
        _buffer = _init_bools(dest, Bool, geom; missingval=false, kw...)
        # Burn each geometry in a loop. This can be threaded because
        # we writ to the same arrays (We could allocate more for that maybe)
        for geom in geoms
            _burn_geometry!(B, geom; _buffer, kw...)
        end
    else
        Threads.@spawn for (i, geom) in enumerate(geoms)
            B1 = view(B, Dim{:geometry}(i))
            _burn_geometry!(B1, geom; kw...)
        end
    end
    return B
end


using PolygonInbounds: vertex, edgecount, edgeindex, searchfirst

# _burn_polygon!
# Fill a raster with `fill` where pixels are inside a polygon
# `boundary` determines how edges are handled
function _burn_polygon!(B::AbstractRaster, geom;
    fill=true, boundary=:center, geomextent, kw...
)
    # Area
    edges, nodes = to_edges_and_nodes(geom)
    B1 = _prepare_for_burning(B, Center())
    _burn_polygon!(B1, nodes, edges)

    # Lines
    n_on_line = false
    if boundary === :touches && _check_intervals(B, boundary)
        # Add line pixels
        B2 = _prepare_for_burning(B, Center())
        n_on_line = _burn_lines!(B2, geom; fill)
    elseif boundary === :inside && _check_intervals(B, boundary)
        # Remove line pixels
        B2 = _prepare_for_burning(B, Center())
        n_on_line = _burn_lines!(B2, geom; fill=!fill)
    elseif boundary !== :center
        throw(ArgumentError("`boundary` can be :touches, :inside, or :center, got $boundary"))
    end
    if verbose
        (n_on_line > 0) || @warn "some polygons were not filled as they did not cross a pixel center. Consider using `boundary=:touches`, or `verbose=false` to hide these warning"
    end
    return B
end

# Modified from PolygonInbounds with optimisations for rasters
# where we know point order and can loop over the array directly
function _burn_polygon!(A::AbstractDimArray, nodes, edges)
    poly = PolygonInbounds.PolygonMesh(nodes, edges)
    nvertices = length(A) # number of points to be checked
    nedges = edgecount(poly) # number of edges of the polygon mesh
    xlookup = lookup(A, XDim)
    ylookup = lookup(A, YDim)
    ix = 1
    iy = 2

    # loop over polygon edges
    for epos = 1:nedges
        inod = edgeindex(poly, epos, 1)  # from
        jnod = edgeindex(poly, epos, 2)  # to
        # swap order of vertices
        if vertex(poly, inod, iy) > vertex(poly, jnod, iy)
            inod, jnod = jnod, inod
        end

        # calc. edge bounding-box
        xone = vertex(poly, inod, ix)
        xtwo = vertex(poly, jnod, ix)
        yone = vertex(poly, inod, iy)
        ytwo = vertex(poly, jnod, iy)

        ymin = yone
        ymax = ytwo

        xmin = min(xone, xtwo)
        xmax = max(xone, xtwo)

        ydel = ytwo - yone
        xdel = xtwo - xone

        # Find Y index 
        ystart = searchsortedfirst(ylookup, ymin)
        # calc. edge-intersection
        # loop over all points with y ∈ [ymin, ymax)
        for y in ystart:LA.ordered_step(ylookup):LA.ordered_lastindex(ylookup)
            @inbounds ypos = ylookup[y]
            ypos > ymax && break 
            for x in eachindex(xlookup)
                @inbounds xpos = xlookup[x]
                xpos > xmax && break 
                if xpos >= xmin
                    # inside extended bounding box of edge
                    mul1 = ydel * (xpos - xone)
                    mul2 = xdel * (ypos - yone)
                    if mul1 < mul2 && yone <= ypos < ytwo
                        # left of line && ypos exact to avoid multiple counting
                        @inbounds A[Y(y), X(x)] = !A[Y(y), X(x)]
                    end
                else # xpos < xmin - left of bounding box
                    if yone <= ypos < ytwo
                        # ypos exact to avoid multiple counting
                        @inbounds A[Y(y), X(x)] = !A[Y(y), X(x)]
                    end
                end
            end
        end
    end
    return A
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

@noinline function _check_intervals(B, boundary)
    if all(map(s -> s isa Intervals, sampling(dims(B)))) 
        return true
    else
        @info "`boundary=:$boundary` only applies to intervals. Line burning skipped"
        return false
    end
end

function _prepare_for_burning(B, locus)
    B1 = _forward_ordered(B)
    start_dims = map(dims(B1, DEFAULT_POINT_ORDER)) do d
        d = DD.maybeshiftlocus(locus, d)
        modify(Array, d)
    end
    return setdims(B1, start_dims)
end

function _forward_ordered(B)
    reduce(dims(B1); init=B1) do A, d
        if DD.order(d) isa ReverseOrdered
            A = view(A, rebuild(d, lastindex(d):-1:firstindex(d)))
            set(A, d => reverse(d))
        else
            A
        end
    end
end


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
    fill, atol=nothing, kw...
)
    selectors = map(dims(x, DEFAULT_POINT_ORDER)) do d
        _at_or_contains(d, _dimcoord(d, point), atol)
    end
    # TODO make a check in dimensionaldata that returns the index if it is inbounds
    if hasselection(x, selectors)
        I = dims2indices(x, selectors)
        _burn_index!(x, fill, I)
    end
    return x
end

# Fill Int indices directly
function _burn_index!(st::AbstractRasterStack, fill, I)
    foreach(st, fill) do A, f
        _burn_index!(A, f, I)
    end
end
_burn_index!(A::AbstractRaster, fill, I::NTuple{<:Any,Int}) = A[I...] = fill
_burn_index!(A::AbstractRaster, fill::Function, I::NTuple{<:Any,Int}) =
    A[I...] = fill(A[I...])
# Handle filling over arbitrary dimensions.
function _burn_index!(A::AbstractRaster, fill, I)
    v = view(A, I...)
    for i in eachindex(v)
        val = fill isa Function ? fill(v[i]) : fill
        v[i] = val
    end
end

# _burn_lines!
# Fill a raster with `fill` where pixels touch lines in a geom
# Separated for a type stability function barrier
function _burn_lines!(B::AbstractRaster, geom; fill=true, kw...)
    # Make sure dims are `Intervals` with `Start` locus
    start_dims = map(dims(B, DEFAULT_POINT_ORDER)) do d
        DD.maybeshiftlocus(Start(), d)
    end
    _prepare_for_burning(B, Start())
    # Flip the order with a view to keep our alg simple
    return _burn_lines!(forward_ordered_B, geom, fill)
end

_burn_lines!(forward_ordered_B, geom, fill) =
    _burn_lines!(forward_ordered_B, GI.geomtrait(geom), geom, fill)
function _burn_lines!(B::AbstractArray, ::Union{GI.MultiLineStringTrait}, geom, fill)
    n_on_line = 0
    for linestring in GI.getlinestring(geom)
        n_on_line += _burn_lines!(B, linestring, fill)
    end
    return n_on_line
end
function _burn_lines!(
    B::AbstractArray, ::Union{GI.MultiPolygonTrait,GI.PolygonTrait}, geom, fill
)
    n_on_line = 0
    for ring in GI.getring(geom)
        n_on_line += _burn_lines!(B, ring, fill)
    end
    return n_on_line
end
function _burn_lines!(
    B::AbstractArray, ::Union{GI.LineStringTrait,GI.LinearRingTrait}, linestring, fill
)
    isfirst = true
    local firstpoint, laststop
    n_on_line = 0
    for point in GI.getpoint(linestring)
        if isfirst
            isfirst = false
            firstpoint = point
            laststop = (x=GI.x(point), y=GI.y(point))
            continue
        end
        if point == firstpoint
            isfirst = true
        end
        line = (
            start=laststop,
            stop=(x=GI.x(point), y=GI.y(point)),
        )
        laststop = line.stop
        n_on_line += _burn_line!(B, line, fill)
    end
    return n_on_line
end
function _burn_lines!(
    B::AbstractArray, t::GI.LineTrait, line, fill
)
    p1, p2 = GI.getpoint(t, line)
    line1 = (
        start=(x=GI.x(p1), y=GI.y(p1)),
        stop=(x=GI.x(p2), y=GI.y(p2)),
    )
    return _burn_line!(B, line1, fill)
end

# _burn_line!
#
# Line-burning algorithm
# Burns a single line into a raster with value where pixels touch a line
#
# TODO: generalise to Irregular spans?
function _burn_line!(A::AbstractRaster, line, fill)
    xdim, ydim = dims(A, DEFAULT_POINT_ORDER)
    regular = map((xdim, ydim)) do d
        lookup(d) isa AbstractSampled && span(d) isa Regular
    end
    msg = """
        Can only fill lines where dimensions have `Regular` lookups.
        Consider using `boundary=:center`, reprojecting the crs,
        or make an issue in Rasters.jl on github if you need this to work.
        """
    all(regular) || throw(ArgumentError(msg))

    @assert order(xdim) == order(ydim) == LookupArrays.ForwardOrdered()
    @assert locus(xdim) == locus(ydim) == LookupArrays.Start()

    raster_x_scale = abs(step(span(A, X)))
    raster_y_scale = abs(step(span(A, Y)))
    raster_x_offset = xdim[1]
    raster_y_offset = ydim[1]
    # Converted lookup to array axis values (still floating)
    relstart = (x=(line.start.x - raster_x_offset) / raster_x_scale, 
             y=(line.start.y - raster_y_offset) / raster_y_scale)
    relstop = (x=(line.stop.x - raster_x_offset) / raster_x_scale, 
            y=(line.stop.y - raster_y_offset) / raster_y_scale)
    diff_x = relstop.x - relstart.x
    diff_y = relstop.y - relstart.y

    # Ray/Slope calculations
    # Straight distance to the first vertical/horizontal grid boundaries
    if relstop.x > relstart.x
        xoffset = floor(relstart.x) - relstart.x + 1 
        xmoves = floor(relstop.x) - floor(relstart.x)
    else
        xoffset = relstart.x - floor(relstart.x)
        xmoves = floor(relstart.x) - floor(relstop.x)
    end
    if relstop.y > relstart.y
        yoffset = floor(relstart.y) - relstart.y + 1
        ymoves = floor(relstop.y) - floor(relstart.y)
    else
        yoffset = relstart.y - floor(relstart.y)
        ymoves = floor(relstart.y) - floor(relstop.y)
    end
    manhattan_distance = xmoves + ymoves
    # Angle of ray/slope.
    # max: How far to move along the ray to cross the first cell boundary.
    # delta: How far to move along the ray to move 1 grid cell.
    hyp = @fastmath sqrt(diff_y^2 + diff_x^2)
    cs = diff_x / hyp
    si = -diff_y / hyp

    delta_x, max_x = if isapprox(cs, zero(cs); atol=1e-10)
        -Inf, Inf
    else
        1.0 / cs, xoffset / cs
    end
    delta_y, max_y = if isapprox(si, zero(si); atol=1e-10)
        -Inf, Inf
    else
        1.0 / si, yoffset / si
    end
    # For arbitrary dimension indexing
    dimconstructors = map(DD.basetypeof, (xdim, ydim))
    # Count how many exactly hit lines
    n_on_line = 0
    countx = county = 0

    # Int starting points for the lin. +1 converts to julia indexing
    j, i = floor(Int, relstart.x) + 1, floor(Int, relstart.y) + 1 # Int

    # Int steps to move allong the line
    step_j = signbit(diff_x) * -2 + 1
    step_i = signbit(diff_y) * -2 + 1

    # Travel one grid cell at a time. Start at zero for the current cell
    for _ in 0:manhattan_distance
        D = map((d, o) -> d(o), dimconstructors, (j, i))
        if checkbounds(Bool, A, D...)
            n_on_line += 1
            val = fill isa Function ? fill(A[D...]) : fill
            @inbounds A[D...] = val
        end

        # Only move in either X or Y coordinates, not both.
        if abs(max_x) < abs(max_y)
            max_x += delta_x
            j += step_j
            countx +=1
        else
            max_y += delta_y
            i += step_i
            county +=1
        end
    end
    return n_on_line
end
function _burn_line!(A::AbstractRaster, line, fill, order::Tuple{Vararg{<:Dimension}})
    msg = """"
        Converting a `:line` geometry to raster is currently only implemented for 2d lines.
        Make a Rasters.jl github issue if you need this for more dimensions.
        """
    throw(ArgumentError(msg))
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
@inline _geom_shape(geom) = _geom_shape(GI.geomtrait(geom), geom)
@inline _geom_shape(::Union{<:GI.PointTrait,<:GI.MultiPointTrait}, geom) = :point
@inline _geom_shape(::Union{<:GI.LineTrait,<:GI.LineStringTrait,<:GI.MultiLineStringTrait}, geom) = :line
@inline _geom_shape(::Union{<:GI.LinearRingTrait,<:GI.PolygonTrait,<:GI.MultiPolygonTrait}, geom) = :polygon
@inline _geom_shape(x, geom) = throw(ArgumentError("Geometry trait $x cannot be rasterized"))
@inline _geom_shape(::Nothing, geom) = throw(ArgumentError("Object is not a GeoInterface.jl compatible geometry: $geom"))

