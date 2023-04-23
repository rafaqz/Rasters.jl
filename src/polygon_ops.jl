const DEFAULT_POINT_ORDER = (X(), Y())
const DEFAULT_TABLE_DIM_KEYS = (:X, :Y)
const INIT_BURNSTATUS = (; ic=1, burn=false, hasburned=false)

# Simple point positions offset relative to the raster
struct Position
    offset::Tuple{Float64,Float64}
    yind::Int32
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
        new(start.offset, gradient, start.yind + 1, stop.yind)
    end
end

function _can_skip_edge(prevpos::Position, nextpos::Position, xlookup, ylookup)
    # ignore horizontal edges
    (prevpos.offset[2] == nextpos.offset[2]) && return true
    # ignore edges between grid lines on y axis
    (nextpos.yind == prevpos.yind) && 
        (prevpos.offset[2] != prevpos.yind) && 
        (nextpos.offset[2] != nextpos.yind) && return true
    # ignore edges outside the grid on the y axis
    (prevpos.offset[2] < 0) && (nextpos.offset[2] < 0) && return true
    (prevpos.offset[2] > lastindex(ylookup)) && (nextpos.offset[2] > lastindex(ylookup)) && return true
    return false
end

struct Allocs{B}
    buffer::B
    edges::Vector{Edge}
    scratch::Vector{Edge}
    crossings::Vector{Float64}
end
Allocs(buffer) = Allocs(buffer, Vector{Edge}(undef, 0), Vector{Edge}(undef, 0), Vector{Float64}(undef, 0))

function _burning_allocs(x; nthreads=_nthreads(), kw...) 
    return [Allocs(_init_bools(x; metadata=Metadata())) for _ in 1:nthreads]
end

_get_alloc(thread_allocs::Vector{<:Allocs}) =
    _get_alloc(thread_allocs[Threads.threadid()])
_get_alloc(allocs::Allocs) = allocs

x_at_y(e::Edge, y) = (y - e.start[2]) * e.gradient + e.start[1]

Base.isless(e1::Edge, e2::Edge) = isless(e1.iystart, e2.iystart)
Base.isless(e::Edge, x::Real) = isless(e.iystart, x)
Base.isless(x::Real, e::Edge) = isless(x, e.iystart)

_to_edges(geom, dims; kw...) = _to_edges(GI.geomtrait(geom), geom, dims; kw...)
function _to_edges(
    tr::Union{GI.LinearRingTrait,GI.AbstractPolygonTrait,GI.AbstractMultiPolygonTrait}, geom, dims;
    allocs::Union{Allocs,Vector{Allocs}}, kw...
)
    (; edges, scratch) = _get_alloc(allocs)
    local edge_count = max_ylen = 0
    if tr isa GI.LinearRingTrait
        edge_count, max_ylen = _to_edges!(edges, geom, dims, edge_count)
    else
        for ring in GI.getring(geom)
             edge_count, ring_max_ylen = _to_edges!(edges, ring, dims, edge_count)
             max_ylen = max(max_ylen, ring_max_ylen)
        end
    end

    # We may have allocated too much
    edges1 = view(edges, 1:edge_count)
    @static if VERSION < v"1.9-alpha1"
        sort!(edges1)
    else
        sort!(edges1; scratch)
    end

    return edges1, max_ylen
end
function _to_edges!(edges, geom, dims, edge_count)
    GI.npoint(geom) > 0 || return edge_count

    # Dummy Initialisation
    local firstpos = prevpos = nextpos = Position((0.0, 0.0), 0)

    firstpoint = true

    local max_ylen = 0

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
        if _can_skip_edge(prevpos, nextpos, xlookup, ylookup)
            prevpos = nextpos
            continue
        end

        # Add the edge to our `edges` vector
        edge_count += 1
        edge = Edge(prevpos, nextpos)
        max_ylen = max(max_ylen, edge.iystop - edge.iystart)
        if edge_count <= lastindex(edges)
            edges[edge_count] = edge
        else
            push!(edges, edge)
        end
        prevpos = nextpos
        prevpoint = p
    end

    return edge_count, max_ylen
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
_burn_geometry!(B, obj; kw...) = _burn_geometry!(B, GI.trait(obj), obj; kw...)::Bool
function _burn_geometry!(B::AbstractRaster, ::GI.AbstractFeatureTrait, feature; kw...)::Bool
    _burn_geometry!(B, GI.geometry(feature); kw...)
end
function _burn_geometry!(B::AbstractRaster, ::GI.AbstractFeatureCollectionTrait, fc; kw...)::Bool
    geoms = (GI.geometry(f) for f in GI.getfeature(gc))
    _burn_geometry!(B, nothing, geoms; kw...)
end
# Where geoms is an iterator
function _burn_geometry!(B::AbstractRaster, trait::Nothing, geoms; 
    collapse::Union{Bool,Nothing}=nothing, lock=SectorLocks(), verbose=true, progress=true, kw...
)::Bool
    thread_allocs = _burning_allocs(B) 
    range = _geomindices(geoms)
    checklock = Threads.SpinLock()
    burnchecks = _alloc_burnchecks(range)
    p = progress ? _progress(length(range)) : nothing
    if isnothing(collapse) || collapse
        Threads.@threads for i in range
            geom = _getgeom(geoms, i)
            ismissing(geom) && continue
            allocs = _get_alloc(thread_allocs)
            B1 = allocs.buffer
            burnchecks[i] = _burn_geometry!(B1, geom; allocs, lock, kw...)
            progress && ProgressMeter.next!(p)
        end
        buffers = map(a -> a.buffer, thread_allocs)
        _do_broadcast!(|, B, buffers...)
    else
        Threads.@threads for i in range
            geom = _getgeom(geoms, i)
            ismissing(geom) && continue
            B1 = view(B, Dim{:geometry}(i))
            allocs = _get_alloc(thread_allocs)
            burnchecks[i] = _burn_geometry!(B1, geom; allocs, lock, kw...)
            progress && ProgressMeter.next!(p)
        end
    end
    
    _set_burnchecks(burnchecks, metadata(B), verbose)
    return false
end

function _burn_geometry!(B::AbstractRaster, ::GI.AbstractGeometryTrait, geom; 
    shape=nothing, verbose=true, boundary=:center, allocs=Allocs(B), kw...
)::Bool
    hasburned = false
    # Use the specified shape or detect it
    shape = shape isa Symbol ? shape : _geom_shape(geom)
    if shape === :point
        hasburned = _fill_point!(B, geom; fill=true, shape, kw...)
    elseif shape === :line
        n_on_line = _burn_lines!(B, geom; shape, kw...)
        hasburned = n_on_line > 0
    elseif shape === :polygon
        # Get the extents of the geometry and array
        geomextent = _extent(geom)
        arrayextent = Extents.extent(B, DEFAULT_POINT_ORDER)
        # Only fill if the gemoetry bounding box overlaps the array bounding box
        if !Extents.intersects(geomextent, arrayextent) 
            verbose && _verbose_extent_info(geomextent, arrayextent)
            return false
        end
        # Take a view of the geometry extent
        B1 = view(B, Touches(geomextent))
        buf1 = _init_bools(B1)
        # Burn the polygon into the buffer
        hasburned = _burn_polygon!(buf1, geom; shape, geomextent, allocs, boundary, kw...)
        for i in eachindex(B1)
            if buf1[i] 
                B1[i] = true
            end
        end
    else
        _shape_error(shape) 
    end
    return hasburned
end

@noinline _shape_error(shape) = 
    throw(ArgumentError("`shape` is $shape, must be `:point`, `:line`, `:polygon` or `nothing`"))

@noinline _verbose_extent_info(geomextent, arrayextent) =
    @info "A geometry was ignored at $geomextent as it was outside of the supplied extent $arrayextent"

# _burn_polygon!
# Burn `true` values into a raster
# `boundary` determines how edges are handled 
function _burn_polygon!(B::AbstractDimArray, geom; kw...)::Bool
    B1 = _prepare_for_burning(B)
    _burn_polygon!(B1::AbstractDimArray, GI.geomtrait(geom), geom; kw...)
end
function _burn_polygon!(B::AbstractDimArray, trait, geom;
    fill=true, boundary=:center, geomextent, verbose=false, allocs=Allocs(B), kw...
)::Bool
    allocs = _get_alloc(allocs)
    filtered_edges, max_ylen = _to_edges(geom, dims(B); allocs)
    
    hasburned::Bool = _burn_polygon!(B, filtered_edges, allocs.crossings, max_ylen)

    # Lines
    n_on_line = 0
    if boundary !== :center
        _check_intervals(B, boundary)
        if boundary === :touches 
            if _check_intervals(B, boundary)
                # Add line pixels
                n_on_line = _burn_lines!(B, geom; fill)::Int
            end
        elseif boundary === :inside 
            if _check_intervals(B, boundary)
                # Remove line pixels
                n_on_line = _burn_lines!(B, geom; fill=!fill)::Int
            end
        else
            throw(ArgumentError("`boundary` can be :touches, :inside, or :center, got :$boundary"))
        end
        if verbose
            (n_on_line > 0) || @info "$n_on_line pixels were on lines"
        end
    end

    hasburned |= (n_on_line > 0)

    return hasburned
end
function _burn_polygon!(A::AbstractDimArray, edges::AbstractArray{Edge}, crossings, max_ylen;
    offset=nothing, verbose=true
)::Bool
    local prev_ypos = 0
    hasburned = false
    # Loop over each index of the y axis
    for iy in axes(A, YDim)
        # Calculate where on the x axis iy is crossed
        ncrossings, prev_ypos = _set_crossings!(crossings, A, edges, iy, prev_ypos, max_ylen)
        # Burn between alternate crossings
        burndata = _burn_crossings!(A, crossings, ncrossings, iy)
        hasburned = burndata.hasburned
    end
    return hasburned
end

function _set_crossings!(crossings, A, edges, iy, prev_ypos, max_ylen)
    ypos = max(1, prev_ypos - max_ylen - 1)
    ncrossings = 0
    # We know the maximum size on y, so we can start from ypos 
    start_ypos = searchsortedfirst(edges, ypos)
    prev_ypos = start_ypos
    for i in start_ypos:lastindex(edges)
        y = iy + first(lookup(A, YDim)) - 1
        e = edges[i]
        # Edges are sorted on y, so we can skip
        # some at the end once they are larger than iy
        if iy < e.iystart 
            prev_ypos = iy
            break
        end
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
    return ncrossings, prev_ypos
end

function _burn_crossings!(A, crossings, ncrossings, iy; 
    status::typeof(INIT_BURNSTATUS)=INIT_BURNSTATUS, 
)
    stop = false
    # Start burning loop from outside any rings
    ic, burn = status
    ix = firstindex(A, X())
    hasburned = false
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
                hasburned = true
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
    return (; ic, burn, hasburned)
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
    return true
end
@noinline function _fill_point!(x::RasterStackOrArray, ::GI.AbstractPointTrait, point;
    fill, atol=nothing, lock=nothing, kw...
)
    selectors = map(dims(x, DEFAULT_POINT_ORDER)) do d
        _at_or_contains(d, _dimcoord(d, point), atol)
    end
    # TODO make a check in dimensionaldata that returns the index if it is inbounds
    if hasselection(x, selectors)
        I = dims2indices(x, selectors)
        if isnothing(lock)  
            _fill_index!(x, fill, I)
        else
            sector = CartesianIndices(map(i -> i:i, I))
            Base.lock(lock, sector)
            _fill_index!(x, fill, I)
            Base.unlock(lock)
        end
        return true
    else
        return false
    end
end

# Fill Int indices directly
_fill_index!(st::AbstractRasterStack, fill::NamedTuple, I::NTuple{<:Any,Int}) = st[I...] = fill
_fill_index!(A::AbstractRaster, fill, I::NTuple{<:Any,Int}) = A[I...] = fill
_fill_index!(A::AbstractRaster, fill::Function, I::NTuple{<:Any,Int}) = A[I...] = fill(A[I...])

_fill_index!(st::AbstractRasterStack, fill::NamedTuple, I) = 
    map((A, f) -> A[I...] .= Ref(f), st, fill)
_fill_index!(A::AbstractRaster, fill, I) = A[I...] .= Ref(fill)
_fill_index!(A::AbstractRaster, fill::Function, I) = A[I...] .= fill.(view(A, I...))

# _burn_lines!
# Fill a raster with `fill` where pixels touch lines in a geom
# Separated for a type stability function barrier
function _burn_lines!(B::AbstractRaster, geom; fill=true, kw...)
    _check_intervals(B)
    B1 = _prepare_for_burning(B)
    return _burn_lines!(B1, geom, fill)
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
    g1 = first(data)
    if GI.trait(g1) isa GI.PointTrait 
        xs = extrema(p -> GI.x(p), data)
        ys = extrema(p -> GI.y(p), data)
        return _float64_xy_extent(Extents.Extent(X=xs, Y=ys))
    else
        ext = reduce(geoms; init=_extent(first(geoms))) do ext, geom
            Extents.union(ext, _extent(geom))
        end
        return _float64_xy_extent(ext)
    end
end
_extent(::Nothing, data::RasterStackOrArray)::XYExtent = _float64_xy_extent(Extents.extent(data))
function _extent(::Nothing, data::T)::XYExtent where T
    if Tables.istable(T)
        geomcolname = first(GI.geometrycolumns(data))
        cols = Tables.columns(data)
        if geomcolname in Tables.columnnames(cols)
            # Table of geometries
            geoms = Tables.getcolumn(data, geomcolname)
            return _extent(nothing, geoms)
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
        ext isa Extents.Extent || throw(ArgumentError("object returns `nothing` from `Extents.extent`."))
        return _float64_xy_extent(ext)
    end
end
function _extent(::GI.AbstractPointTrait, point)::XYExtent
    x, y = Float64(GI.x(point)), Float64(GI.y(point))
    Extents.Extent(X=(x, x), Y=(y, y))
end
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
_init_bools(to; kw...) = _init_bools(to, Bool; kw...)
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
function _init_bools(to, dims::DimTuple, T::Type, data; collapse::Union{Bool,Nothing}=nothing, kw...)
    if isnothing(data) || isnothing(collapse) || collapse
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

function _alloc_bools(to, dims::DimTuple, ::Type{Bool}; missingval=false, metadata=NoMetadata(), kw...)
    # Use a BitArray
    return Raster(falses(size(dims)), dims; missingval, metadata) # Use a BitArray
end
function _alloc_bools(to, dims::DimTuple, ::Type{T}; missingval=false, metadata=NoMetadata(), kw...) where T
    # Use an `Array`
    data = fill!(Raster{T}(undef, dims), missingval) 
    return rebuild(data; missingval, metadata)
end

_alloc_burnchecks(n::Int) = fill(false, n)
_alloc_burnchecks(x::AbstractArray) = _alloc_burnchecks(length(x))
function _set_burnchecks(burnchecks, metadata::Metadata{<:Any,<:Dict}, verbose)
    metadata["missed_geometries"] = .!burnchecks
    verbose && _burncheck_info(burnchecks)
end
_set_burnchecks(burnchecks, metadata, verbose) = verbose && _burncheck_info(burnchecks)
function _burncheck_info(burnchecks)
    nburned = sum(burnchecks)
    nmissed = length(burnchecks) - nburned
    nmissed > 0 && @info "$nmissed geometries did not affect any pixels. See `metadata(raster)[\"missed_geometries\"]` for a vector of misses"
end

_nthreads() = Threads.nthreads()

