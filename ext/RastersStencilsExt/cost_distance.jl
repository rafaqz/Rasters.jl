function _initialise_cost!(active, accumulated_costs, origins, hood)
    for I in CartesianIndices(origins)
        o = @inbounds origins[I]
        if !ismissing(o) && o < Inf * u"hr"
            @inbounds accumulated_costs[I] = o
            # Add an index to the `active` set if it has any neighbors not
            # in `origins` as we only need to use the edges of the origin areas.
            if any(x -> !ismissing(x) && x == Inf * u"hr", Stencils.stencil(origins, I))
                push!(active, I)
            end
        end
    end
    return active
end

"""
    cost_distance([f=walking_speed]; origins, elevation, [missingval])

Calculate the cost-distance between cells in `origins` and
all other points in the `origins`/`costs`.

`origins` is a matrix of `Real` where values larger than `1` are used
as starting points with cost of zero.

`costs` is a matrix where values to pass to `a` and `b` in function `f`

The function `f` of the form `f(elevation1, elevation2, distance, resistance)`

- `a`: The value in `costs` at the current location.
- `b`: The value in `costs` at the next location.
- `distance`: The distance between `a` and `b` locations.

"""
cost_distance(origins::AbstractRaster, elevation::AbstractRaster, resistance=nothing; kw...)
    cost_distance(time_taken, origins, elevation, resistance; kw...)
function cost_distance(
    f::Function, origins::AbstractRaster, elevation::AbstractRaster, resistance=nothing; 
    kw...
)
    # Create a grid we will update cost distance into.
    # The initial cost is the maximum `Float64` value. This will later
    # be the `missingval`, as it will remain in any cells that were not touched by
    # the algorithm (because one of the other arrays has missing values there).
    stencil = Moore{1}()
    accumulated_time = StencilArray(fill(Inf * u"hr", size(origins)), stencil; padding=Conditional(), boundary=Remove(typemax(Float64)))
    origins = StencilArray(parent(origins), stencil; padding=Conditional(), boundary=Remove(typemax(Float64)))
    # Write the cost-distance to the accumulator grid
    cost_distance!(accumulated_time, origins, elevation, resistance; stencil, kw...)
    # Remove the padding edge, and wrap the output as an AbstractRaster
    # if costs was one, updating the name and `missingval`.
    if elevation isa AbstractRaster
        return rebuild(elevation; data=parent(accumulated_time), name=:time, missingval=Inf * u"hr")
    else
        return Stencils.unpad_view(acc, 1)
    end
end
function cost_distance!(f, accumulated_time, origins, elevation, resistance=nothing; 
    cellsize=1, missingval=missingval(elevation), stencil
)
    # The neighborood is a simple 3 * 3 moore neighborhood (ring)
    # The active cells are a `Set` of `CartesianIndices` linking
    # to positions in acc/origins/costs.
    # "active" means assigned their lowest cost in the last round,
    # but with neighbors whose cost has not been recalulated.
    active_cells = Set{CartesianIndex{2}}()
    new_active_cells = Set{CartesianIndex{2}}()
    # We add the first active cells from the `origins` array
    _initialise_cost!(active_cells, accumulated_time, origins, stencil)
    # And calculate how many active cells there are
    n_active_cells = length(active_cells)
    # Now loop while there are still active cells that are
    # reducing cost-distance somewhere on the grid.
    while n_active_cells > 0
        for I in active_cells
            # Get the current accumulated cost
            @inbounds current_time = accumulated_time[I]
            # Missing cells are skipped
            e1 = elevation[I]
            ismissingval(e1, missingval) && continue
            # Loop over the neighborhood offsets and distances from center cell
            for (O, d) in zip(Neighborhoods.cartesian_offsets(stencil), Neighborhoods.distances(stencil))
                # Get the index of the neighboring cell
                NI = O + I
                checkbounds(Bool, elevation, NI) || continue
                e2 = elevation[NI]
                # Out of bounds cells are skipped because our costs are not padded
                # Missing valued neighbors are skipped
                ismissing(e2) && continue
                @inbounds current_neighbor_time = accumulated_time[NI]
                # Calculate the new time by adding the time to get to the
                # neighbor to the time to get to the current cell
                r = isnothing(resistance) ? 1 : resistance[I]
                comb = f(e1, e2, d * cellsize, r)
                # @show comb
                new_neighbor_time = current_time + comb
                # @show new_neighbor_time current_neighbor_time
                # Update the time if the new time is lower than any previous path
                if new_neighbor_time < current_neighbor_time
                    @inbounds accumulated_time[NI] = new_neighbor_time
                    push!(new_active_cells, NI)
                end
            end
        end
        # Remove the active cells, their neighbors have been calculated and
        # added to new_active_cells, so they are done unless a cheaper path reaches
        # them again later on.
        empty!(active_cells)
        # Swap the active and new active Dicts for the next round of calculations
        active_cells, new_active_cells = new_active_cells, active_cells
        # Update how many cells we are calculating next round
        n_active_cells = length(active_cells)
    end
    return accumulated_time
end

ismissingval(val, missingval) = val === missingval
ismissingval(vals::NamedTuple, missingval) = any(map(val -> ismissingval(val, missingval), vals))
ismissingval(vals::NamedTuple, missingvals::NamedTuple) = any(map(ismissingval, vals, missingvals))

"""
    walking_time(e1, e2, d, v)

Calculate the cost of moving between elevation `e1` to elevation `e2`
over distance `d` through resistance `r`

cost = S * I * A * E
"""
function walking_time(e1, e2, d, r)
    s = walking_speed(ImhofTobler(), abs(e1 - e2), d)
    return d / (s * r)
end

abstract type AnisotropicCost end

struct ImhofTobler <: AnisotropicCost
    slopefactor::Float64
    distfactor::Float64
end
ImhofTobler(; slopefactor=-3.5, distfactor=0.05) =
    ImhofTobler(slopefactor, distfactor)

walking_speed(x::ImhofTobler, dvert, dhoriz) = walking_speed(x, dvert / dhoriz)
walking_speed(x::ImhofTobler, slope) = (6ℯ^(x.slopefactor * (slope + x.distfactor))) * u"km/hr"

abstract type IsotropicCost end

struct SouleGoldman <: IsotropicCost end

land_cover_speed(::SouleGoldman, category) = relative_movement_speed[category]

const movement_speed = (
    road=5.0u"km/h",
    dirt_road=4.0u"km/h",
    track=3.2u"km/h",
    cleared_land=2.0u"km/h",
    forest=1.0u"km/h",
    dense_forest=0.5u"km/h",
    wetland=0.5u"km/h",
)

const relative_movement_speed = map(movement_speed) do s
    s / movement_speed.road
end


