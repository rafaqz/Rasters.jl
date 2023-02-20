const DEFAULT_POINT_ORDER = (X(), Y())
const DEFAULT_TABLE_DIM_KEYS = (:X, :Y)

# Simple point positions offset relative to the raster
struct Position
    offset::Tuple{Float64,Float64}
    ind::Int32
end
function Position(point::Tuple, start::Tuple, step::Tuple)
    (x, y) = point
    xoff, yoff = (x - start[1]), (y - start[2])
    offset = (xoff / step[1] + 1.0, yoff / step[2] + 1.0)
    yind = trunc(Int32, offset[2])
    return Position(offset, yind)
end

# Simple edge offset with (x, y) start relative to the raster
# gradient of line and integer start/stop for y
struct Edge
    start::Tuple{Float64,Float64}
    gradient::Float64
    iystart::Int32
    iystop::Int32
    function Edge(start::Position, stop::Position)
        if start.offset[2] > stop.offset[2]
            stop, start = start, stop
        end
        gradient = (stop.offset[1] - start.offset[1]) / (stop.offset[2] - start.offset[2])
        new(start.offset, gradient, start.ind + 1, stop.ind)
    end
end

function can_skip_edge(prevpos::Position, nextpos::Position, xlookup, ylookup)
    # ignore horizontal edges
    (prevpos.offset[2] == nextpos.offset[2]) && return true
    # ignore edges between grid lines on y axis
    (nextpos.ind == prevpos.ind) && 
        (prevpos.offset[2] != prevpos.ind) && 
        (nextpos.offset[2] != nextpos.ind) && return true
    # ignore edges outside the grid on the y axis
    (prevpos.offset[2] < 0) && (nextpos.offset[2] < 0) && return true
    (prevpos.offset[2] > lastindex(ylookup)) && (nextpos.offset[2] > lastindex(ylookup)) && return true
    return false
end

struct Allocs
    edges::Vector{Edge}
    scratch::Vector{Edge}
    crossings::Vector{Float64}
    spinlock::Threads.SpinLock
end
Allocs() = Allocs(0)
Allocs(n::Int) = Allocs(Vector{Edge}(undef, n), Vector{Edge}(undef, n), Vector{Float64}(undef, n), Threads.SpinLock())

_edge_allocs() = _edge_allocs(Threads.nthreads())
_edge_allocs(n::Int) = [Allocs() for _ in 1:n]

_get_alloc(thread_allocs::Vector{Allocs}; kw...) =
    _get_alloc(thread_allocs[Threads.threadid()])
_get_alloc(allocs::Allocs; kw...) = allocs

x_at_y(e::Edge, y) = (y - e.start[2]) * e.gradient + e.start[1]

Base.isless(e1::Edge, e2::Edge) = isless(e1.iystart, e2.iystart)


_to_edges(geom, dims; kw...) = _to_edges(GI.geomtrait(geom), geom, dims; kw...)
function _to_edges(
    tr::Union{GI.LinearRingTrait,GI.AbstractPolygonTrait,GI.AbstractMultiPolygonTrait}, geom, dims;
    allocs::Union{Allocs,Vector{Allocs}}, kw...
)
    (; edges, scratch) = _get_alloc(allocs)
    local edge_count = 0
    if tr isa GI.LinearRingTrait
        edge_count = _to_edges!(edges, geom, dims, edge_count)
    else
        for ring in GI.getring(geom)
             edge_count = _to_edges!(edges, ring, dims, edge_count)
        end
    end

    # We may have allocated too much
    edges1 = view(edges, 1:edge_count)
    sort!(edges1; scratch)
    return edges1
end
function _to_edges!(edges, geom, dims, edge_count)
    GI.npoint(geom) > 0 || return edge_count

    # Dummy Initialisation
    local firstpos = prevpos = nextpos = Position((0.0, 0.0), 0)

    firstpoint = true

    # Raster properties
    xlookup, ylookup = lookup(dims, (X(), Y())) 
    starts = (Float64(first(xlookup)), Float64(first(ylookup)))
    steps = (Float64(Base.step(xlookup)), Float64(Base.step(ylookup)))
    local prevpoint = (0.0, 0.0)

    # Loop over points to generate edges
    for point in GI.getpoint(geom)
        p = (Float64(GI.x(point)), Float64(GI.y(point)))
       
        # For the first point just set variables
        if firstpoint
            prevpos = firstpos = Position(p, starts, steps)
            prevpoint = p
            firstpoint = false
            continue
        end

        # Get the next offsets and indices
        nextpos = Position(p, starts, steps)

        # Check if we need an edge between these offsets
        # This is the performance-critical step that reduces the size of the edge list
        if can_skip_edge(prevpos, nextpos, xlookup, ylookup)
            prevpos = nextpos
            continue
        end

        # Add the edge to our `edges` vector
        edge_count += 1
        edge = Edge(prevpos, nextpos)
        if edge_count <= lastindex(edges)
            edges[edge_count] = edge
        else
            push!(edges, edge)
        end
        prevpos = nextpos
        prevpoint = p
    end

    return edge_count
end

# _burn_geometry!
# Fill a raster with `fill` where it interacts with a geometry.
# This is used in `boolmask` TODO move to mask.jl ?
# 
# _istable keyword is a hack so we know not to pay the
# price of calling `istable` which calls `hasmethod`
function burn_geometry!(B::AbstractRaster, data::T; kw...) where T
    if Tables.istable(T)
        geomcolname = first(GI.geometrycolumns(data))::Symbol
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
    shape=nothing, _buffer=nothing, verbose=false, allocs, kw...
)
    allocs = _get_alloc(allocs)
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
            verbose && _verbose_extent_info(geomextent, arrayextent)
            return B
        end
        # Take a view of the geometry extent
        B1 = view(B, Touches(geomextent))
        buf1 = if isnothing(_buffer)
            _init_bools(B1, Bool; kw..., missingval=false)
        else
            buf1 = view(_buffer, Touches(geomextent))
        end
        # Also take a view of the buffer
        # Reset the buffer for this area
        fill!(buf1, false)
        # Burn the polygon into the buffer
        _burn_polygon!(buf1, geom; shape, geomextent, allocs, kw...)

        # We are writing to the same array with all threads
        lock(allocs.spinlock)
        for i in eachindex(B1)
            if buf1[i] 
                B1[i] = true
            end
        end
        unlock(allocs.spinlock)
    else
        _shape_error(shape) 
    end
    return B
end

@noinline _shape_error(shape) = 
    throw(ArgumentError("`shape` is $shape, must be `:point`, `:line`, `:polygon` or `nothing`"))

@noinline _verbose_extent_info(geomextent, arrayextent) =
    @info "A geometry was ignored at $geomextent as it was outside of the supplied extent $arrayextent"
# Treat geoms as an iterator
function _burn_geometry!(B::AbstractRaster, trait::Nothing, geoms; combine::Union{Bool,Nothing}=nothing, kw...)
    allocs = _edge_allocs()
    if isnothing(combine) || combine
        p = _progress(length(geoms))
        Threads.@threads :static for i in _geomindices(geoms)
            geom = _getgeom(geoms, i)
            ismissing(geom) && continue
            _burn_geometry!(B, geom; allocs, kw...)
            ProgressMeter.next!(p)
        end
    else
        p = _progress(length(geoms))
        Threads.@threads :static for i in _geomindices(geoms)
            geom = _getgeom(geoms, i)
            ismissing(geom) && continue
            B1 = view(B, Dim{:geometry}(i))
            _burn_geometry!(B1, geom; kw...)
            ProgressMeter.next!(p)
        end
    end
    return B
end


# _burn_polygon!
# Burn `true` values into a raster
# `boundary` determines how edges are handled 
function _burn_polygon!(B::AbstractDimArray, geom; kw...)
    B1 = _prepare_for_burning(B)
    _burn_polygon!(B1::AbstractDimArray, GI.geomtrait(geom), geom; kw...)
end
function _burn_polygon!(B::AbstractDimArray, trait, geom;
    fill=true, boundary=:center, geomextent, verbose=false, allocs, kw...
)
    allocs = _get_alloc(allocs)
    filtered_edges = _to_edges(geom, dims(B); allocs)
    
    _burn_polygon!(B, filtered_edges, allocs.crossings)

    # Lines
    n_on_line = false
    if boundary !== :center
        _check_intervals(B, boundary)
        if boundary === :touches && _check_intervals(B, boundary)
            # Add line pixels
            n_on_line = _burn_lines!(B, geom; fill)
        elseif boundary === :inside && _check_intervals(B, boundary)
            # Remove line pixels
            n_on_line = _burn_lines!(B, geom; fill=!fill)
        else
            throw(ArgumentError("`boundary` can be :touches, :inside, or :center, got :$boundary"))
        end
    end
    if verbose
        (n_on_line > 0) || @warn "some polygons were not filled as they did not cross a pixel center. Consider using `boundary=:touches`, or `verbose=false` to hide these warning"
    end
    return B
end

function _burn_polygon!(A::AbstractDimArray, edges::AbstractArray{Edge}, crossings;
    offset= nothing
)
    # Loop over each index of the y axis
    for iy in axes(A, YDim)
        ncrossings = _set_crossings!(crossings, A, edges, iy)
        _burn_crossings!(A, crossings, ncrossings, iy)
    end
    return A
end

function _set_crossings!(crossings, A, edges, iy)
    ncrossings = 0
    # Edges are sorted on y, so we can skip
    # some at the start we have already done
    for i in 1:lastindex(edges)
        y = iy + first(lookup(A, YDim)) - 1
        e = edges[i]
        iy >= e.iystart || break
        if iy <= e.iystop 
            ncrossings += 1
            if length(crossings) >= ncrossings
                crossings[ncrossings] = Rasters.x_at_y(e, iy)
            else
                push!(crossings, Rasters.x_at_y(e, iy))
            end
        end
    end
    sort!(view(crossings, 1:ncrossings), )
    return ncrossings
end

function _burn_crossings!(A, crossings, ncrossings, iy; 
    status::Tuple{Int,Bool}=(1, false), 
)
    stop = false
    # Start burning loop from outside any rings
    ic, burn = status
    ix = firstindex(A, X())
    while ic <= ncrossings
        crossing = crossings[ic]
        # Burn/skip until we hit the next edge crossing
        while ix < crossing
            if ix > lastindex(A, X()) 
                stop = true
                break
            end
            if burn
                A[X(ix), Y(iy)] = true
            end
            ix += 1
        end
        if stop
            break
        else
            # Alternate burning/skipping with each edge crossing
            burn = !burn
            ic += 1
        end
    end
    # Maybe fill in the end of the row
    if burn
        for x in ix:lastindex(A, X())
            A[X(ix), Y(iy)] = true
        end
    end
    return ic, burn
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
        _lookup_as_array(d)
    end
    return setdims(B1, start_dims)
end

# Convert to Array if its not one already
_lookup_as_array(d::Dimension) = parent(lookup(d)) isa Array ? d : modify(Array, d) 

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

    # TODO merge this with Edge generation
    # Converted lookup to array axis values (still floating)
    relstart = (x=(line.start.x - raster_x_offset) / raster_x_step, 
             y=(line.start.y - raster_y_offset) / raster_y_step)
    relstop = (x=(line.stop.x - raster_x_offset) / raster_x_step, 
            y=(line.stop.y - raster_y_offset) / raster_y_step)

    # Ray/Slope calculations
    # Straight distance to the first vertical/horizontal grid boundaries
    if relstop.x > relstart.x
        xoffset = trunc(relstart.x) - relstart.x + 1 
        xmoves = trunc(Int, relstop.x) - trunc(Int, relstart.x)
    else
        xoffset = relstart.x - trunc(relstart.x)
        xmoves = trunc(Int, relstart.x) - trunc(Int, relstop.x)
    end
    if relstop.y > relstart.y
        yoffset = trunc(relstart.y) - relstart.y + 1
        ymoves = trunc(Int, relstop.y) - trunc(Int, relstart.y)
    else
        yoffset = relstart.y - trunc(relstart.y)
        ymoves = trunc(Int, relstart.y) - trunc(Int, relstop.y)
    end
    manhattan_distance = xmoves + ymoves

    # Int starting points for the line. +1 converts to julia indexing
    j, i = trunc(Int, relstart.x) + 1, trunc(Int, relstart.y) + 1 # Int

    # For arbitrary dimension indexing
    dimconstructors = map(DD.basetypeof, (xdim, ydim))

    if manhattan_distance == 0
        D = map((d, o) -> d(o), dimconstructors, (j, i))
        if checkbounds(Bool, A, D...)
            @inbounds A[D...] = fill
        end
        n_on_line = 1
        return n_on_line
    end

    diff_x = relstop.x - relstart.x
    diff_y = relstop.y - relstart.y

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
    # Count how many exactly hit lines
    n_on_line = 0
    countx = county = 0


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

const XYExtent = Extents.Extent{(:X,:Y),Tuple{Tuple{Float64,Float64},Tuple{Float64,Float64}}}

# _extent
# Get the bounds of a geometry
_extent(geom)::XYExtent = _extent(GI.trait(geom), geom)
function _extent(::Nothing, data::AbstractVector)::XYExtent
    ext = reduce(data; init=_extent(first(data))) do ext, geom
        Extents.union(ext, _extent(geom))
    end
    return _float64_xy_extent(ext)
end
_extent(::Nothing, data::RasterStackOrArray)::XYExtent = _float64_xy_extent(Extents.extent(data))
function _extent(::Nothing, data::T)::XYExtent where T
    if Tables.istable(T)
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
            bounds = reduce(DEFAULT_TABLE_DIM_KEYS; init=(;)) do acc, key
                if key in Tables.columnnames(cols)
                    merge(acc, (; key=extrema(cols[key])))
                else
                    acc
                end
            end
            return _float64_xy_extent(Extent(bounds))
        end
    else
        ext = Extents.extent(data)
        ext isa Extent || throw(ArgumentError("object returns `nothing` from `Extents.extent`."))
        return _float64_xy_extent(ext)
    end
end
_extent(::GI.AbstractPointTrait, geom)::XYExtent = Extents.Extent(X=GI.x(geom), Y=GI.y(geom))
function _extent(::GI.AbstractGeometryTrait, geom)::XYExtent
    geomextent = GI.extent(geom)
    if isnothing(geomextent)
        points = GI.getpoint(geom)
        xbounds = extrema(GI.x(p) for p in points)
        ybounds = extrema(GI.y(p) for p in points)
        return _float64_xy_extent(Extents.Extent(X=xbounds, Y=ybounds))
    else
        return _float64_xy_extent(geomextent)
    end
end
_extent(::GI.AbstractFeatureTrait, feature)::XYExtent = _extent(GI.geometry(feature))
function _extent(::GI.AbstractFeatureCollectionTrait, features)::XYExtent
    features = GI.getfeature(features)
    init = _float64_xy_extent(_extent(first(features)))
    ext = reduce(features; init) do acc, f
        Extents.union(acc, _extent(f))
    end
    return _float64_xy_extent(ext)
end

function _float64_xy_extent(ext::Extents.Extent)
    xbounds = map(Float64, ext.X)
    ybounds = map(Float64, ext.Y)
    return Extents.Extent(X=xbounds, Y=ybounds)
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


# Like `create` but without disk writes, mostly for Bool/Union{Missing,Boo},
# and uses `similar` where possible
# TODO merge this with `create` somehow
_init_bools(to, T::Type; kw...) = _init_bools(to, T, nothing; kw...)
_init_bools(to::AbstractRasterSeries, T::Type, data; kw...) = _init_bools(first(to), T, data; kw...)
_init_bools(to::AbstractRasterStack, T::Type, data; kw...) = _init_bools(first(to), T, data; kw...)
_init_bools(to::AbstractRaster, T::Type, data; kw...) = _init_bools(to, dims(to), T, data; kw...)
_init_bools(to::Extents.Extent, T::Type, data; kw...) = _init_bools(to, _extent2dims(to; kw...), T, data; kw...)
_init_bools(to::DimTuple, T::Type, data; kw...) = _init_bools(to, to, T, data; kw...)
function _init_bools(to::Nothing, T::Type, data; kw...)
    # Get the extent of the geometries
    ext = _extent(data)
    isnothing(ext) && throw(ArgumentError("no recognised dimensions, extent or geometry"))
    # Convert the extent to dims (there must be `res` or `size` in `kw`)
    dims = _extent2dims(ext; kw...)
    return _init_bools(to, dims, T, data; kw...)
end
function _init_bools(to, dims::DimTuple, T::Type, data; combine::Union{Bool,Nothing}=nothing, kw...)
    if isnothing(data) || isnothing(combine) || combine
        _alloc_bools(to, dims, T; kw...)
    else
        n = if Base.IteratorSize(data) isa Base.HasShape
            length(data)
        else
            count(_ -> true, data)
        end
        geomdim = Dim{:geometry}(1:n)
        _alloc_bools(to, (dims..., geomdim), T; kw...)
    end
end

function _alloc_bools(to, dims::DimTuple, ::Type{Bool}; missingval, kw...)
    # Use a BitArray
    return Raster(falses(size(dims)), dims; missingval) # Use a BitArray
end
function _alloc_bools(to, dims::DimTuple, ::Type{T}; missingval, kw...) where T
    # Use an `Array`
    data = fill!(Raster{T}(undef, dims), missingval) 
    return Raster(data, dims; missingval)
end
