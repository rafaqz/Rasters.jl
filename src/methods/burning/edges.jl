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

# Edge offset with (x, y) start relative to the raster
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

Base.isless(e1::Edge, e2::Edge) = isless(e1.iystart, e2.iystart)
Base.isless(e::Edge, x::Real) = isless(e.iystart, x)
Base.isless(x::Real, e::Edge) = isless(x, e.iystart)

_x_at_y(e::Edge, y) = (y - e.start[2]) * e.gradient + e.start[1]

#= Edges: a ollection of Edge
- `max_ylen` tracks the longest y-y distance in any edge.
    we can use this to reduce the search space later.
- `edge_count` counts how many edges we have - it may be 
    less than the length of the vector
=#
struct Edges <: AbstractVector{Edge}
    edges::Vector{Edge}
    max_ylen::Int
    min_y::Int
    edge_count::Int
end
Edges(geom, dims; kw...) = Edges(GI.geomtrait(geom), geom, dims; kw...)
function Edges(
    tr::Union{GI.AbstractCurveTrait,GI.AbstractPolygonTrait,GI.AbstractMultiPolygonTrait}, 
    geom, dims;
    allocs, kw...
)
    (; edges, scratch) = _get_alloc(allocs)

    # TODO fix bug that requires this to be redefined
    edges = Vector{Edge}(undef, 0)
    local edge_count = max_ylen = 0
    local min_y = typemax(Int)
    if tr isa GI.AbstractCurveTrait
        edge_count, max_ylen, min_y = _to_edges!(edges, geom, dims, edge_count)
    else
        for ring in GI.getring(geom)
             edge_count, ring_max_ylen, ring_min_y = _to_edges!(edges, ring, dims, edge_count)
             max_ylen = max(max_ylen, ring_max_ylen)
             min_y = min(min_y, ring_min_y)
        end
    end

    # We may have allocated too much
    edges1 = view(edges, 1:edge_count)
    @static if VERSION < v"1.9-alpha1"
        sort!(edges1)
    else
        sort!(edges1; scratch)
    end

    return Edges(edges, max_ylen, min_y, edge_count)
end

Base.parent(edges::Edges) = edges.edges
Base.length(edges::Edges) = edges.edge_count
Base.axes(edges::Edges) = axes(parent(edges))
Base.getindex(edges::Edges, I...) = getindex(parent(edges), I...)
Base.setindex(edges::Edges, x, I...) = setindex!(parent(edges), x, I...)

function _can_skip(prevpos::Position, nextpos::Position, xlookup, ylookup)
    # ignore edges between grid lines on y axis
    (nextpos.yind == prevpos.yind) && 
        (prevpos.offset[2] != prevpos.yind) && 
        (nextpos.offset[2] != nextpos.yind) && return true
    # ignore edges outside the grid on the y axis
    (prevpos.offset[2] < 0) && (nextpos.offset[2] < 0) && return true
    (prevpos.offset[2] > lastindex(ylookup)) && (nextpos.offset[2] > lastindex(ylookup)) && return true
    # ignore horizontal edges
    (prevpos.offset[2] == nextpos.offset[2]) && return true
    return false
end

@noinline function _to_edges!(edges, geom, dims, edge_count)
    # Dummy Initialisation
    local firstpos = prevpos = nextpos = Position((0.0, 0.0), 0)
    isfirst = true
    local max_ylen = 0
    local min_y = typemax(Int)

    GI.npoint(geom) > 0 || return edge_count, max_ylen, min_y
    xlookup, ylookup = lookup(dims, (X(), Y())) 
    (length(xlookup) > 0 && length(ylookup) > 0) || return edge_count, max_ylen, min_y

    # Raster properties
    starts = (Float64(first(xlookup)), Float64(first(ylookup)))
    steps = (Float64(Base.step(xlookup)), Float64(Base.step(ylookup)))
    local prevpoint = (0.0, 0.0)

    # Loop over points to generate edges
    for point in GI.getpoint(geom)
        p = (Float64(GI.x(point)), Float64(GI.y(point)))
       
        # For the first point just set variables
        if isfirst
            prevpos = firstpos = Position(p, starts, steps)
            prevpoint = p
            isfirst = false
            continue
        end

        # Get the next offsets and indices
        nextpos = Position(p, starts, steps)

        # Check if we need an edge between these offsets
        # This is the performance-critical step that reduces the size of the edge list
        if _can_skip(prevpos, nextpos, xlookup, ylookup)
            prevpos = nextpos
            continue
        end

        # Add the edge to our `edges` vector
        edge_count += 1
        edge = Edge(prevpos, nextpos)
        _add_edge!(edges, edge, edge_count)
        max_ylen = max(max_ylen, edge.iystop - edge.iystart)
        min_y = min(min_y, edge.iystart)
        prevpos = nextpos
        prevpoint = p
    end
    # Check in case the polygon is not closed
    if prevpos != firstpos
        edge_count += 1
        edge = Edge(prevpos, firstpos)
        # Update the longest y distance of any edge
        max_ylen = max(max_ylen, edge.iystop - edge.iystart)
        min_y = min(min_y, edge.iystart)
        # assign/push the edge to edges
        _add_edge!(edges, edge, edge_count)
    end

    return edge_count, max_ylen, min_y
end

function _add_edge!(edges, edge, edge_count)
    if edge_count <= lastindex(edges)
        @inbounds edges[edge_count] = edge
    else
        push!(edges, edge)
    end
end
