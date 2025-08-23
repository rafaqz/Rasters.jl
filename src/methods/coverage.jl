const COVERAGE_DOC = """
Calculate the area of a raster covered by GeoInterface.jl compatible geometry `geom`,
as a fraction.

Each pixel is assigned a grid of points (by default 10 x 10) that are each checked
to be inside the geometry. The sum divided by the number of points to give coverage.

In practice, most pixel coverage is not calculated this way - shortcuts that 
produce the same result are taken wherever possible.

If `geom` is an `AbstractVector` or table, the `mode` keyword will determine how coverage is combined.
"""

const COVERAGE_KEYWORDS = """
- `mode`: method for combining multiple geometries - `union` or `sum`. 
    * `union` (the default) gives the areas covered by all geometries. Usefull in
      spatial coverage where overlapping regions should not be counted twice.
      The returned raster will contain `Float64` values between `0.0` and `1.0`.
    * `sum` gives the summed total of the areas covered by all geometries,
      as in taking the sum of running `coverage` separately on all geometries.
      The returned values are positive `Float64`.
    For a single geometry, the `mode` keyword has no effect - the result is the same.
- `scale`: `Integer` scale of pixel subdivision. The default of `10` means each pixel has 
    10 x 10 or 100 points that contribute to coverage. Using `100` means 10,000 points
    contribute. Performance will decline as `scale` increases. Memory use will grow 
    by `scale^2` when `mode=:union`.
$THREADED_KEYWORD
$PROGRESS_KEYWORD
$VERBOSE_KEYWORD
"""

"""
    coverage(mode, geom; [to, res, size, scale, verbose, progress])
    coverage(geom; [to, mode, res, size, scale, verbose, progress])

$COVERAGE_DOC

# Keywords

$COVERAGE_KEYWORDS
$TO_KEYWORD
$GEOMETRYCOLUMN_KEYWORD
$SIZE_KEYWORD
$RES_KEYWORD
"""
coverage(data; to=nothing, mode=union, scale=10,geometrycolumn=nothing, kw...) = _coverage(to, data; mode, scale, geometrycolumn, kw...)
coverage(f::Union{typeof(sum),typeof(union)}, data; geometrycolumn=nothing, kw...) = coverage(data; geometrycolumn, kw..., mode=f)

function _coverage(to, data; mode, scale, kw...)
    name = if GI.isgeometry(data) || GI.isfeature(data)
        Symbol(:coverage)
    else
        Symbol(:coverage_, mode)
    end
    # We use sum for `reducer` so eltype inference works
    r = Rasterizer(data; reducer=sum, fill=0.0, init=0.0, missingval=0.0, kw...)
    rc = RasterCreator(to, data; kw..., eltype=r.eltype, fill, name, missingval=r.missingval)
    dest = create_rasterize_dest(rc) do A
        _coverage!(A, r; mode, scale)
    end
    return dest
end

"""
    coverage!(A, geom; [mode, scale])

$COVERAGE_DOC

# Keywords

$COVERAGE_KEYWORDS
"""
coverage!(mode::Union{typeof(union),typeof(sum)}, A::AbstractRaster, data; kw...) =
    coverage!(A, data; kw..., mode)
function coverage!(A::AbstractRaster, data; scale::Integer=10, mode=union, kw...)
    # We use sum for `reducer` so eltype inference works
    r = Rasterizer(data; reducer=sum, fill=0.0, init=0.0, missingval=0.0, kw...)
    coverage!(A, r; scale, mode) 
end
# Collect iterators so threading is easier.
function _coverage!(A::AbstractRaster, r::Rasterizer; scale::Integer=10, mode=union)  
    _coverage!(A, GI.trait(r.geom), r.geom, r; scale, mode)
end
function _coverage!(A::AbstractRaster, ::GI.AbstractGeometryTrait, geom, r; scale, mode)
    subpixel_dims = _subpixel_dims(A, scale)
    missed_pixels = if mode === union
        _union_coverage!(A, geom; scale, subpixel_dims, progress=r.progress)
    elseif mode === sum
        _sum_coverage!(A, geom; scale, subpixel_dims, progress=r.progress)
    else
        throw(ArgumentError("Coverage `mode` can be `union` or `sum`. Got $mode"))
    end
    r.verbose && _check_missed_pixels(missed_pixels, scale)
    return A
end
function _coverage!(A::AbstractRaster, ::Nothing, geoms, r::Rasterizer; mode, scale)
    n = r.threaded ? _nthreads() : 1
    buffers = (
        # We always need a vector of threads, but just
        allocs = _burning_allocs(A; nthreads=n, threaded=true),
        linebuffer = [_init_bools(A, BitArray; missingval=false) for _ in 1:n],
        centerbuffer = [_init_bools(A, BitArray; missingval=false) for _ in 1:n],
        block_crossings = [[Vector{Float64}(undef, 0) for _ in 1:scale] for _ in 1:n],
        burnstatus=[fill(BurnStatus(), scale) for _ in 1:n],
        subbuffer = [fill!(Array{Bool}(undef, scale, scale), false) for _ in 1:n],
        ncrossings = [fill(0, scale) for _ in 1:n],
    )
    subpixel_dims = _subpixel_dims(A, scale)
    missed_pixels = if mode === union
        _union_coverage!(A, geoms, buffers; scale, subpixel_dims, progress=r.progress, threaded=r.threaded)
    elseif mode === sum
        _sum_coverage!(A, geoms, buffers; scale, subpixel_dims, progress=r.progress, threaded=r.threaded)
    else
        throw(ArgumentError("Coverage `mode` can be `union` or `sum`. Got $mode"))
    end
    r.verbose && _check_missed_pixels(missed_pixels, scale)
    return A
end

# Combines coverage at the sub-pixel level for a final value 0-1
function _union_coverage!(A::AbstractRaster, geoms, buffers; 
    scale, subpixel_dims, progress=true, threaded=false, checkmem=CHECKMEM[],
)
    threaded, n = _check_buffer_thread_mem(A; scale, threaded, checkmem)
    buffer_size = size(A) .* scale
    centeracc = [_init_bools(A, BitArray; missingval=false) for _ in 1:n]
    lineacc = [_init_bools(A, BitArray; missingval=false) for _ in 1:n]
    subpixel_buffer = [falses(buffer_size) for _ in 1:n]
    allbuffers = merge(buffers, (; centeracc, lineacc, subpixel_buffer))
    range = _geomindices(geoms)
    _run(range, threaded, progress, "Calculating coverage buffers...") do i
        geom = _getgeom(geoms, i)
        idx = Threads.threadid()
        # Get buffers for each thread as a NamedTuple
        thread_buffers = map(b -> b[idx], allbuffers)
        _union_coverage!(A, geom; scale, subpixel_buffer, thread_buffers...)
        fill!(thread_buffers.linebuffer, false)
        fill!(thread_buffers.centerbuffer, false) # Is this necessary?
    end

    # Merge downscaled BitArray (with a function barrier)
    subpixel_union = _do_broadcast!(|, subpixel_buffer[1], subpixel_buffer...)
    subpixel_raster = Raster(subpixel_union, subpixel_dims)
    # Merge main BitArray (with a function barrier)
    center_covered = _do_broadcast!(|, centeracc[1], centeracc...)
    line_covered = _do_broadcast!(|, lineacc[1], lineacc...)

    missed_pixels = fill(0, n)
    range = axes(A, Y())
    _run(range, threaded, progress, "Combining coverage...") do y
        for x in axes(A, X())
            D = (X(x), Y(y))
            if center_covered[D...]
                pixel_coverage = 1.0
                A[D...] = pixel_coverage
            else
                subxs, subys = map((x, y)) do i
                    i1 = (i - 1) * scale
                    sub_indices = i1 + 1:i1 + scale
                end
                pixel_coverage = sum(view(subpixel_raster, X(subxs), Y(subys))) / scale^2
                if pixel_coverage == 0.0
                    # Check if this line should be covered
                    if line_covered[D...]
                        idx = Threads.threadid()
                        missed_pixels[idx] += 1
                    end
                end
                A[D...] = pixel_coverage
            end
        end
    end
    return sum(missed_pixels)
end
function _union_coverage!(A::AbstractRaster, geom;
    scale,
    linebuffer=_init_bools(A, BitArray; missingval=false),
    centerbuffer=_init_bools(A, BitArray; missingval=false),
    allocs=Allocs(linebuffer),
    subpixel_buffer=nothing,
    centeracc=_init_bools(A, BitArray; missingval=false),
    lineacc=_init_bools(A, BitArray; missingval=false),
    burnstatus=[BurnStatus() for _ in 1:scale],
    block_crossings=[Vector{Float64}(undef, 0) for _ in 1:scale],
    subbuffer=falses(scale, scale),
    subpixel_dims=_subpixel_dims(A, scale),
    ncrossings=fill(0, scale),
    checkmem=CHECKMEM[],
)

    if isnothing(subpixel_buffer)
        checkmem && _check_buffer_mem(A, scale)
        subpixel_buffer = falses(size(A) .* scale)
    end
    GI.isgeometry(geom) || error("not a geometry")
    crossings = allocs.crossings
    boolmask!(linebuffer, geom; shape=:line, allocs)
    boolmask!(centerbuffer, geom; boundary=:center, allocs)
    # Update the cumulative state for completely covered cells
    # This allows us to skip them later.
    # We don't just use `boundary=:inside` because we need the line burns separately anyway
    centeracc .|= (centerbuffer .& .!(linebuffer))
    lineacc .|= linebuffer
    edges = Edges(geom, subpixel_dims; allocs)

    prev_ypos = 0
    # Loop over y in A
    for y in axes(A, Y())
        # If no lines touched this column skip it
        found = false
        for x in axes(A, X())
            if linebuffer[X(x), Y(y)] && !centeracc[X(x), Y(y)]
                found = true
                break
            end
        end
        found || continue

        y1 = (y - 1) * scale
        sub_yaxis = y1 + 1:y1 + scale

        # Generate all of the x crossings beforehand so we don't do it for every pixel
        for (i, sub_y) in enumerate(sub_yaxis)
            ncrossings[i], prev_ypos = _set_crossings!(block_crossings[i], edges, sub_y, prev_ypos)
        end

        # Reset burn burnstatus
        for i in eachindex(burnstatus)
            burnstatus[i] = BurnStatus()
        end
        subpixel_raster = Raster(subpixel_buffer, subpixel_dims)

        # Loop over x in A
        for x in axes(A, X())
            # Checkt that the cell is not already 100% filled
            centeracc[X(x), Y(y)] && continue
            # CHeck that a line touched the cell
            linebuffer[X(x), Y(y)] || continue
            x1 = (x - 1) * scale
            sub_xaxis = x1 + 1:x1 + scale
            offset_subbuffer = OffsetArrays.OffsetArray(subbuffer, (sub_xaxis, sub_yaxis))
            subdims = map(subpixel_dims, axes(offset_subbuffer)) do d, a
                rebuild(d, NoLookup(a)) # Don't need a lookup for _burn_polygon
            end

            # Rebuild the buffer for this pixels dims
            subraster = Raster(offset_subbuffer, subdims, (), NoName(), NoMetadata(), false)
            # And initialise it
            fill!(subraster, false)
            # Loop over y in the subraster
            for (i, sub_y) in enumerate(sub_yaxis)
                ncrossings[i] > 0 || continue
                # Burn along x for each y, tracking burn status and crossing number
                burnstatus[i] = _burn_crossings!(subraster, block_crossings[i], ncrossings[i], sub_y; status=burnstatus[i])
            end
            v = view(subpixel_raster, X(sub_xaxis), Y(sub_yaxis))
            parent(v) .|= parent(parent(subraster))
        end
    end
end

# Sums all coverage 
function _sum_coverage!(A::AbstractRaster, geoms, buffers; 
    scale, subpixel_dims, verbose=true, progress=true, threaded=false
)
    n = _nthreads()
    coveragebuffers = [fill!(similar(A), 0.0) for _ in 1:n]
    missed_pixels = fill(0, n)
    range = _geomindices(geoms)
    burnchecks = _alloc_burnchecks(range)
    _run(range, threaded, progress, "Calculating coverage...") do i
        geom = _getgeom(geoms, i)
        ismissing(geom) && return nothing
        idx = Threads.threadid()
        thread_buffers = map(b -> b[idx], buffers)
        coveragebuffer = coveragebuffers[idx]
        nmissed, burnchecks[i] = _sum_coverage!(coveragebuffer, geom; scale, thread_buffers...)
        missed_pixels[idx] += nmissed
        fill!(thread_buffers.linebuffer, false)
        fill!(thread_buffers.centerbuffer, false)
        return nothing
    end
    _do_broadcast!(+, A, coveragebuffers...)
    _set_burnchecks(burnchecks, metadata(A), verbose)
    return sum(missed_pixels)
end
function _sum_coverage!(A::AbstractRaster, geom;
    scale,
    linebuffer=_init_bools(A, BitArray; missingval=false),
    centerbuffer=_init_bools(A, BitArray; missingval=false),
    allocs=Allocs(linebuffer),
    block_crossings=[Vector{Float64}(undef, 0) for _ in 1:scale],
    subbuffer=falses(scale, scale),
    burnstatus=[BurnStatus() for _ in 1:scale],
    subpixel_dims=_subpixel_dims(A, scale),
    ncrossings=fill(0, scale),
)
    GI.isgeometry(geom) || error("Object is not a geometry")
    crossings = allocs.crossings
    boolmask!(linebuffer, geom; shape=:line, allocs)
    boolmask!(centerbuffer, geom; boundary=:center, allocs)
    edges = Edges(geom, subpixel_dims; allocs)
    # Brodcast over the rasterizations and indices
    # to calculate coverage of each pixel
    local missed_pixels = 0
    local hascover = false

    prev_ypos = 0
    # Loop over y in A
    for y in axes(A, Y())
        found = false
        for x in axes(A, X())
            if linebuffer[X(x), Y(y)]
                found = true
            elseif centerbuffer[X(x), Y(y)]
                # If the center is inside a polygon but the pixel is
                # not on a line, then coverage is 1.0
                pixel_coverage = 1.0
                hascover = true
                A[X(x), Y(y)] += pixel_coverage
            end
        end
        # If no lines touched this column, skip it
        found || continue

        y1 = (y - 1) * scale
        sub_yaxis = y1 + 1:y1 + scale

        # Generate all of the x crossings beforehand so we don't do it for every pixel
        for (i, sub_y) in enumerate(sub_yaxis)
            ncrossings[i], prev_ypos = _set_crossings!(block_crossings[i], edges, sub_y, prev_ypos)
        end
        # Set the burn/skip status to false (skip) for each starting position
        burnstatus .= Ref(BurnStatus())

        missed_pixels = 0
        # Loop over x in A
        for x in axes(A, X())
            linebuffer[X(x), Y(y)] || continue
            hasburned =

            x1 = (x - 1) * scale
            sub_xaxis = x1 + 1:x1 + scale
            offset_subbuffer = OffsetArrays.OffsetArray(subbuffer, (sub_xaxis, sub_yaxis))
            subdims = map(subpixel_dims, axes(offset_subbuffer)) do d, a
                # Don't need a real lookup for _burn_polygon
                rebuild(d, NoLookup(a)) 
            end

            # Rebuild the buffer for this pixels dims
            subraster = Raster(offset_subbuffer, subdims, (), NoName(), NoMetadata(), false)
            # And initialise it
            fill!(subraster, false)

            # Loop over y in the subraster
            for (i, sub_y) in enumerate(sub_yaxis)
                ncrossings[i] > 0 || continue
                # Burn along x for each y, tracking burn status and crossing number
                burnstatus[i] = _burn_crossings!(subraster, block_crossings[i], ncrossings[i], sub_y; status=burnstatus[i])
            end

            # Finally, sum the pixel sub-raster and divide by scale^2 for the fractional coverage
            pixel_coverage = sum(subraster) / scale^2
            if pixel_coverage == 0
                missed_pixels += 1
            else
                hascover = true
            end
            A[X(x), Y(y)] += pixel_coverage
        end
    end
    return missed_pixels, hascover 
end

function _subpixel_dims(A, scale)
    shifted = map(d -> DD.maybeshiftlocus(Start(), d), commondims(A, DEFAULT_TABLE_DIM_KEYS))
    map(shifted) do d
        l = lookup(d)
        substep = step(l) / scale
        substart = substep / 2 + first(l)
        range = substart:substep:last(l) + scale * substep
        sublookup = Sampled(range, ForwardOrdered(), Regular(substep), Intervals(Start()), NoMetadata())
        rebuild(d, sublookup)
    end
end

_buffer_bytes(A, scale) = prod(size(A) .* scale) ÷ 8 # bits to bytes

function _check_buffer_thread_mem(A; scale, threaded, checkmem)
    n = threaded ? _nthreads() : 1
    if checkmem && n > 1 && Sys.free_memory() < _buffer_bytes(A, scale) * n
        @warn "Not enough memory to use `threaded=true` with `scale=$scale`. Using `threaded=false`"
        threaded = false
        n = 1
    end
    checkmem && _check_buffer_mem(A, scale)
    return threaded, n
end

function _check_buffer_mem(A, scale)
    bytes = _buffer_bytes(A, scale)
    f = x -> let scale=scale
        "Not enough memory for `coverage` at `scale=$scale`. Try a smaller number for the `scale` keyword."
    end
    _checkmem(f, bytes)
end

function _check_missed_pixels(missed_pixels::Int, scale::Int)
    if missed_pixels > 0
        @info "There were $missed_pixels times that pixels were touched by geometries but did not produce coverage, meaning these areas did not contain any of the $(scale * scale) points in the pixel at `scale=$scale`."
    end
end
