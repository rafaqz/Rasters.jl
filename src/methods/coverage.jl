"""
    coverage(geoms; to, kw...)
   
Calculate the area of a raster covered by geometries, as a fraction.

If `geoms` is an `AbstractVector` or table, coverage will be summed,
giving results more than 1.
"""
coverage(data; to=nothing, kw...) = _coverage(to, data; kw...)

_coverage(to::Extents.Extent, data; res=nothing, size=nothing, kw...) =
    _coverage(_extent2dims(to; res, size, kw...), data; kw...)
_coverage(to::Nothing, data; kw...) = _coverage(_extent(data), data; kw...)
function _coverage(to, data; kw...)
    dest = _create_rasterize_dest(dims(to); fill=0.0, init=0.0, name=:coverage, missingval=0.0, kw...) do dest
        coverage!(dest, data; kw...)
    end
    return dest
end

coverage!(A::AbstractRaster, data; kw...) = coverage!(A, GI.trait(data), data; kw...)
function coverage!(dest, ::Nothing, data::AbstractVector; kw...)
    n = Threads.nthreads()
    thread_allocs = _edge_allocs() 
    linebuffers = [_init_bools(dest, Bool; missingval=false) for _ in 1:n]
    centerbuffers = [_init_bools(dest, Bool; missingval=false) for _ in 1:n]
    coveragebuffers = [fill!(similar(dest), 0.0) for _ in 1:n]
    subbuffers = [fill!(Array{Bool}(undef, 10, 10), false) for _ in 1:n]
    Threads.@threads for geom in data
        idx = Threads.threadid()
        coveragebuffer = coveragebuffers[idx]
        linebuffer = linebuffers[idx]
        centerbuffer = centerbuffers[idx]
        subbuffer = subbuffers[idx]
        allocs = thread_allocs[idx]
        coverage!(coveragebuffer, geom; allocs, linebuffer, centerbuffer, subbuffer)
        fill!(linebuffer, false)
        fill!(centerbuffer, false)
    end
    dest .= .+(coveragebuffers...)
    return dest
end
function coverage!(A::AbstractRaster, ::Union{GI.PolygonTrait,GI.MultiPolygonTrait}, data;
    allocs=Allocs(),
    linebuffer=_init_bools(A, Bool; missingval=false),
    centerbuffer=_init_bools(A, Bool; missingval=false),
    subbuffer=falses(10, 10),
)
    crossings = allocs.crossings
    lines = boolmask!(linebuffer, data; shape=:line, allocs)
    centers = boolmask!(centerbuffer, data; boundary=:center, allocs)
    shifted = map(d -> DD.maybeshiftlocus(Start(), d), dims(A, DEFAULT_TABLE_DIM_KEYS))
    disag_dims = map(shifted) do d
        l = lookup(d)
        substep = step(l) / 10
        substart = substep / 2 + first(l)
        range = substart:substep:last(l) + 10 * substep
        sublookup = Sampled(range, ForwardOrdered(), Regular(substep), Intervals(Start()), NoMetadata())
        rebuild(d, sublookup)
    end
    filtered_edges = _to_edges(data, disag_dims; allocs) |> sort!


    # Brodcase over the rasterizations and indices
    # to calculate coverage of each pixel
    edgestart = 1
    block_crossings = [Vector{Float64}(undef, 0) for _ in 1:10]
    status = [(1, false) for _ in 1:10]

    # Loop over y in A
    for y in axes(A, Y())
        # If the center is in a polygon, and the pixel is
        # not on a line, then coverage is 1.0
        for x in axes(A, X())
            if centers[X(x), Y(y)] && !lines[X(x), Y(y)] 
                A[X(x), Y(y)] += 1.0
            end
        end

        # If no lines touched this column skip it
        any(view(lines, Y(y))) || continue

        y1 = (y - 1) * 10
        sub_yaxis = y1 + 1:y1 + 10

        # Generate all of the x crossings beforehand so we don't do it for every pixel
        for (i, sub_y) in enumerate(sub_yaxis)
            edgestart, ncrossings = _set_crossings!(crossings, A, filtered_edges, edgestart, sub_y)
            block_crossings[i] = crossings[1:ncrossings] # Allocation :(
        end
        status .= Ref((1, false)) 

        # Loop over x in A
        for x in axes(A, X())
            lines[X(x), Y(y)] || continue
            x1 = (x - 1) * 10
            sub_xaxis = x1 + 1:x1 + 10
            offset_subbuffer = OffsetArrays.OffsetArray(subbuffer, (sub_xaxis, sub_yaxis))
            subdims = map(disag_dims, axes(offset_subbuffer)) do d, a
                rebuild(d, NoLookup(a)) # Don't need a lookup for _burn_polygon
            end

            # Rebuild the buffer for this pixels dims
            subraster = Raster(offset_subbuffer, subdims, (), NoName(), NoMetadata(), false)
            # And initialise it
            fill!(subraster, false)
            # Loop over y in the subraster
            for (i, sub_y) in enumerate(sub_yaxis)
                cross = block_crossings[i]
                ncrossings = length(cross)
                ncrossings > 0 || continue
                # Burn along x for each y, tracking burn status and crossing number
                status[i] = _burn_crossings!(subraster, cross, ncrossings, sub_y; status=status[i])
            end
            # Sum the raster and divide by 100 for fractional coverage
            pixel_coverage = sum(subraster) / 100.0
            A[X(x), Y(y)] += pixel_coverage
        end
    end
end
