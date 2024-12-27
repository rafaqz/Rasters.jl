struct Allocs{B}
    buffer::B
    edges::Vector{Edge}
    scratch::Vector{Edge}
    crossings::Vector{Float64}
end
function Allocs(buffer)
    edges = Vector{Edge}(undef, 0)
    scratch = Vector{Edge}(undef, 0)
    crossings = Vector{Float64}(undef, 0)
    return Allocs(buffer, edges, scratch, crossings)
end

function _burning_allocs(x::Nothing; 
    nthreads=_nthreads(), 
    threaded=true, 
    burncheck_metadata=nothing,
    kw...
) 
    threaded ? [Allocs(nothing) for _ in 1:nthreads] : [Allocs(nothing)]
end
function _burning_allocs(x; 
    nthreads=_nthreads(), 
    threaded=true, 
    burncheck_metadata=Metadata(),
    kw...
) 
    dims = commondims(x, (XDim, YDim))
    if threaded
        [Allocs(_init_bools(dims; metadata=deepcopy(burncheck_metadata))) for _ in 1:nthreads]
    else
        [Allocs(_init_bools(dims; metadata=burncheck_metadata))]
    end
end

function _get_alloc(allocs::Vector{<:Allocs}) 
    if length(allocs) == 1
        only(allocs)
    else
        allocs[Threads.threadid()]
    end
end
_get_alloc(allocs::Allocs) = allocs

# TODO include these in Allocs
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
