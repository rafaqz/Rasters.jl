
const DEFAULT_POINT_ORDER = (X(), Y())
const DEFAULT_TABLE_DIM_KEYS = (:X, :Y)

# Moved here form PolygonInbounds for optimisation
# Thanks to Klaus Crusius for the code here and in _burn_polygon
struct PolygonMesh{U,E}
    nodes::U
    edges::E
end
PolygonMesh(geom) = PolygonMesh(_to_nodes_and_edges(geom)...)

nodecount(poly::PolygonMesh) = size(poly.nodes, 1)
Base.length(poly::PolygonMesh) = nodecount(poly)
edgecount(poly::PolygonMesh) = length(poly.edges)
edgeindex(poly::PolygonMesh, i::Integer, n::Integer) = poly.edges[i][n]

vertex(poly::PolygonMesh, v::Integer, xy::Integer) = poly.nodes[v][xy]

_to_nodes_and_edges(geom) = _to_nodes_and_edges(GI.geomtrait(geom), geom)
function _to_nodes_and_edges(
    tr::Union{GI.LinearRingTrait,GI.AbstractPolygonTrait,GI.AbstractMultiPolygonTrait}, geom
)
    n = GI.npoint(geom)
    edges = Vector{Tuple{Int,Int}}(undef, n)
    nodes = Vector{Tuple{Float64,Float64}}(undef, n)
    nodenum = 0
    if tr isa GI.LinearRingTrait
        _to_edges_and_nodes!(edges, nodes, nodenum, geom)
    else
        for ring in GI.getring(geom)
            nodenum = _to_nodes_and_edges!(edges, nodes, nodenum, ring)
        end
    end
    return nodes, edges 
end
function _to_nodes_and_edges!(edges, nodes, lastnode, geom)
    npoints = GI.npoint(geom)
    for (n, point) in enumerate(GI.getpoint(geom))
        i = lastnode + n
        if n == npoints
            # The closing edge of a sub-polygon
            edges[i] = (i, lastnode + 1)
        else
            # A regular edge somewhere in a sub-polygon
            edges[i] = (i, i + 1)
        end
        nodes[i] = (GI.x(point), GI.y(point))
    end
    return lastnode + npoints
end


struct Edge{Float64}
    start::Tuple{T,T}
    stop::Tuple{T,T}
    istart::Tuple{Int,Int}
    istop::Tuple{Int,Int}
    closing::Bool
end

_to_edges(geom, dims) = _to_edges(GI.geomtrait(geom), geom, dims)
function _to_edges(
    tr::Union{GI.LinearRingTrait,GI.AbstractPolygonTrait,GI.AbstractMultiPolygonTrait}, geom, dims
)
    edges = Vector{Tuple{Int,Int}}(undef, GI.npoint(geom))
    if tr isa GI.LinearRingTrait
        _to_edges!(edges, geom)
    else
        for ring in GI.getring(geom)
             _to_edges!(edges, ring)
        end
    end
    return edges 
end
function _to_edges!(edges, geom, dims)
    xlookup, ylookup = lookup(dims, (XDim, YDim)) 
    xstart, ystart = first(xlookup), first(ylookup)
    xstep, ystep = step(xlookup), step(ylookup)

    npoints = GI.npoint(geom)
    first = true
    for (n, point) in enumerate(GI.getpoint(geom))
        if first
            firstpoint = lastpoint = x, y = (GI.x(point), GI.y(point))
            firstind = lastind = floor(Int, (x - xstart) / xstep), floor(Int, (y - ystart) / ystep) 
            first = false
            continue
        end
        # A regular edge somewhere in a sub-polygon
        nextpoint = x, y = (GI.x(point), GI.y(point))
        nextind = floor(Int, (x - xstart) / xstep), floor(Int, (y - ystart) / ystep) 
        edges[i - 1] = Edge(lastpoint, nextpoint, lastind, nextind, false)
        if n == npoints
            edges[i] = Edge(nextpoint, firstpoint, nextind, firstind, false)
            firstpoint = x, y = (GI.x(point), GI.y(point))
            nextind = floor(Int, (x - xstart) / xstep), floor(Int, (y - ystart) / ystep) 
        end
        lastpoint = nextpoint
        lastind = nextind
    end
    return nothing
end


# _burn_geometry!
# Fill a raster with `fill` where it interacts with a geometry.
# This is used in `boolmask` TODO move to mask.jl ?
# 
# _istable keyword is a hack so we know not to pay the
# price of calling `istable` which calls `hasmethod`
function burn_geometry!(B::AbstractRaster, data; kw...)
    if Tables.istable(typeof(data))
        geomcolname = first(GI.geometrycolumns(data))
        geoms = Tables.getcolumn(data, geomcolname)
        _burn_geometry!(B, nothing, geoms; kw...)
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
    geoms = (GI.geometry(f) for f in GI.getfeature(gc))
    _burn_geometry!(B, nothing, geoms; kw...)
end
function _burn_geometry!(B::AbstractRaster, ::GI.AbstractGeometryTrait, geom; 
    shape=nothing, _buffer=nothing, verbose=false, kw...
)
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
            buf1 = _init_bools(B1, Bool, nothing; kw..., missingval=false)
            _burn_polygon!(buf1, geom; shape, geomextent, kw...)
            for i in eachindex(B1)
                if buf1[i] 
                    B1[i] = true # Only writing true is immune to race conditions
                end
            end
        else
            # Also take a view of the buffer
            buf1 = view(_buffer, Touches(geomextent))
            # Reset the buffer for this area
            buf1 .= false
            # Burn the polygon into the buffer
            _burn_polygon!(buf1, geom; shape, geomextent, kw...)
            for i in eachindex(B1)
                if buf1[i] 
                    B1[i] = true # Only writing true is immune to race conditions
                end
            end
        end
    else
        throw(ArgumentError("`shape` is $shape, must be `:point`, `:line`, `:polygon` or `nothing`"))
    end
    return B
end
# Treat geoms as an iterator
function _burn_geometry!(B::AbstractRaster, trait::Nothing, geoms; combine=true, kw...)
    task = if combine
        Threads.@threads for geom in geoms
            _burn_geometry!(B, geom; kw...)
        end
    else
        Threads.@threads for (i, geom) in enumerate(geoms)
            B1 = view(B, Dim{:geometry}(i))
            _burn_geometry!(B1, geom; kw...)
        end
    end
    return B
end


# _burn_polygon!
# Fill a raster with `fill` where pixels are inside a polygon
# `boundary` determines how edges are handled
function _burn_polygon!(B::AbstractDimArray, geom;
    fill=true, boundary=:center, geomextent, verbose=false, kw...
)
    # Area
    B1 = _prepare_for_burning(B)
    mesh = PolygonMesh(geom)
    _burn_polygon!(B1, mesh)

    # Lines
    n_on_line = false
    if boundary !== :center
        _check_intervals(B, boundary)
        if boundary === :touches && _check_intervals(B, boundary)
            # Add line pixels
            n_on_line = _burn_lines!(B1, geom; fill)
        elseif boundary === :inside && _check_intervals(B, boundary)
            # Remove line pixels
            n_on_line = _burn_lines!(B1, geom; fill=!fill)
        else
            throw(ArgumentError("`boundary` can be :touches, :inside, or :center, got :$boundary"))
        end
    end
    if verbose
        (n_on_line > 0) || @warn "some polygons were not filled as they did not cross a pixel center. Consider using `boundary=:touches`, or `verbose=false` to hide these warning"
    end
    return B
end

# Modified from PolygonInbounds with optimisations for rasters
# where we know point order and can loop over the array directly
function _burn_polygon!(A::AbstractDimArray, mesh::PolygonMesh)
    nvertices = length(A) # number of points to be checked
    nedges = edgecount(mesh) # number of edges of the polygon mesh
    xlookup = lookup(A, XDim)
    ylookup = lookup(A, YDim)
    ix = 1
    iy = 2
    sort(

    # loop over polygon edges
    for epos = 1:nedges
        inod = edgeindex(mesh, epos, 1)  # from
        jnod = edgeindex(mesh, epos, 2)  # to
        # swap order of vertices
        if vertex(mesh, inod, iy) > vertex(mesh, jnod, iy)
            inod, jnod = jnod, inod
        end

        # calc. edge bounding-box
        xone = vertex(mesh, inod, ix)
        xtwo = vertex(mesh, jnod, ix)
        yone = vertex(mesh, inod, iy)
        ytwo = vertex(mesh, jnod, iy)

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
        for y in ystart:lastindex(ylookup)
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

const INTERVALS_INFO = "makes more sense on `Intervals` than `Points` and will have more correct results. You can construct dimensions with a `X(values; sampling=Intervals(Center()))` to acheive this"

@noinline _check_intervals(B) = 
    _chki(B) ? true : (@info "burning lines $INTERVALS_INFO"; false)
@noinline _check_intervals(B, boundary) =
    _chki(B) ? true : (@info "`boundary=:$boundary` $INTERVALS_INFO"; false)

_chki(B) = all(map(s -> s isa Intervals, sampling(dims(B)))) 

function _prepare_for_burning(B, locus=Center())
    B1 = _forward_ordered(B)
    start_dims = map(dims(B1, DEFAULT_POINT_ORDER)) do d
        # Shift lookup values to center of pixels
        d = DD.maybeshiftlocus(locus, d)
        # Convert to Array if its not one already
        parent(lookup(d)) isa Array ? d : modify(Array, d) 
    end
    return setdims(B1, start_dims)
end

function _forward_ordered(B)
    reduce(dims(B); init=B) do A, d
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
_burn_index!(st::AbstractRasterStack, fill::Union{Tuple,NamedTuple}, I) = st[I...] = fill
_burn_index!(A::AbstractRaster, fill, I::NTuple{<:Any,Int}) = A[I...] = fill
_burn_index!(A::AbstractRaster, fill::Function, I::NTuple{<:Any,Int}) =
    A[I...] = fill(A[I...])
# # Handle filling over arbitrary dimensions.
# function _burn_index!(A::AbstractRaster, fill, I)
#     v = view(A, I...)
#     for i in eachindex(v)
#         val = fill isa Function ? fill(v[i]) : fill
#         v[i] = val
#     end
# end

# _burn_lines!
# Fill a raster with `fill` where pixels touch lines in a geom
# Separated for a type stability function barrier
function _burn_lines!(B::AbstractRaster, geom; fill=true, kw...)
    _check_intervals(B)
    B1 = _prepare_for_burning(B)
    _burn_lines!(B1, geom, fill)
    return B
end

_burn_lines!(B, geom, fill) =
    _burn_lines!(B, GI.geomtrait(geom), geom, fill)
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
        @assert (parent(lookup(d)) isa Array)
        lookup(d) isa AbstractSampled && span(d) isa Regular
    end
    msg = """
        Can only fill lines where dimensions have `Regular` lookups.
        Consider using `boundary=:center`, reprojecting the crs,
        or make an issue in Rasters.jl on github if you need this to work.
        """
    all(regular) || throw(ArgumentError(msg))

    @assert order(xdim) == order(ydim) == LookupArrays.ForwardOrdered()
    @assert locus(xdim) == locus(ydim) == LookupArrays.Center()

    raster_x_step = abs(step(span(A, X)))
    raster_y_step = abs(step(span(A, Y)))
    raster_x_offset = @inbounds xdim[1] - raster_x_step / 2 # Shift from center to start of pixel
    raster_y_offset = @inbounds ydim[1] - raster_y_step / 2
    # Converted lookup to array axis values (still floating)
    relstart = (x=(line.start.x - raster_x_offset) / raster_x_step, 
             y=(line.start.y - raster_y_offset) / raster_y_step)
    relstop = (x=(line.stop.x - raster_x_offset) / raster_x_step, 
            y=(line.stop.y - raster_y_offset) / raster_y_step)
    diff_x = relstop.x - relstart.x
    diff_y = relstop.y - relstart.y

    # Ray/Slope calculations
    # Straight distance to the first vertical/horizontal grid boundaries
    if relstop.x > relstart.x
        xoffset = floor(relstart.x) - relstart.x + 1 
        xmoves = floor(Int, relstop.x) - floor(Int, relstart.x)
    else
        xoffset = relstart.x - floor(relstart.x)
        xmoves = floor(Int, relstart.x) - floor(Int, relstop.x)
    end
    if relstop.y > relstart.y
        yoffset = floor(relstart.y) - relstart.y + 1
        ymoves = floor(Int, relstop.y) - floor(Int, relstart.y)
    else
        yoffset = relstart.y - floor(relstart.y)
        ymoves = floor(Int, relstart.y) - floor(Int, relstop.y)
    end
    manhattan_distance = xmoves + ymoves
    # Angle of ray/slope.
    # max: How far to move along the ray to cross the first cell boundary.
    # delta: How far to move along the ray to move 1 grid cell.
    hyp = @fastmath sqrt(diff_y^2 + diff_x^2)
    cs = diff_x / hyp
    si = -diff_y / hyp

    delta_x, max_x =# if isapprox(cs, zero(cs); atol=1e-10)
        # -Inf, Inf
    # else
        1.0 / cs, xoffset / cs
    # end
    delta_y, max_y =# if isapprox(si, zero(si); atol=1e-10)
        # -Inf, Inf
    # else
        1.0 / si, yoffset / si
    # end
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
            @inbounds A[D...] = fill
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
    geomextent = GI.extent(geom; )
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

