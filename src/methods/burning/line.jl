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
    B::AbstractArray, ::GI.AbstractCurveTrait, linestring, fill
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
