
"""
    reproject(obj; crs)

Reproject the lookups of `obj` to a different crs. 

This is a lossless operation for the raster data, as only the 
lookup values change. This is only possible when the axes of source
and destination projections are alligned: the change is usually from
a [`Regular`](@ref) and an [`Irregular`](@ref) lookup spans.

For converting between projections that are rotated, 
skewed or warped in any way, `resample` use [`resample`](@ref).

Dimensions without an `AbstractProjected` lookup (such as a `Ti` dimension)
are silently returned without modification.

# Arguments

- `obj`: a `LookupArray`, `Dimension`, `Tuple` of `Dimension`, `Raster` or `RasterStack`.
$CRS_KEYWORD
"""
reproject(x; crs::GeoFormat) = reproject(crs, x)
reproject(target::GeoFormat, x) = rebuild(x; dims=reproject(target, dims(x)))
reproject(target::GeoFormat, dims::Tuple) = map(d -> reproject(target, d), dims)
reproject(target::GeoFormat, l::LookupArray) = l
reproject(target::GeoFormat, dim::Dimension) = rebuild(dim, reproject(target, lookup(dim)))
function reproject(target::GeoFormat, l::AbstractProjected)
    source = crs(l)
    newdata = reproject(source, target, l.dim, parent(l))
    newlookup = rebuild(l; data=newdata, crs=target)
    if isregular(newdata)
        return set(newlookup, Regular(stepof(newdata)))
    else
        newbounds = reproject(crs(l), target, l.dim, bounds(l))
        return set(newlookup, Irregular(newbounds))
    end
end

"""
    reproject(source::GeoFormat, target::GeoFormat, dim::Dimension, val)

`reproject` uses ArchGDAL.reproject, but implemented for a reprojecting
a value array of values, a single dimension at a time.
"""
function reproject(source, target, dim, val)
    if source == target
        return val
    else
        return _reproject(source, target, dim, val)
    end
end

_reproject(source::GeoFormat, target::GeoFormat, dim, val) = val
_reproject(source::Nothing, target::GeoFormat, dim, val) = val
_reproject(source::GeoFormat, target::Nothing, dim, val) = val
_reproject(source::Nothing, target::Nothing, dim, val) = val

function _reproject(source::GeoFormat, target::GeoFormat, dim::Union{XDim,YDim}, val::Number)
    # This is a dumb way to do this. But it save having to inspect crs, 
    # and prevents reprojections that dont make sense from working.
    # A better method for this should be implemented in future.
    return first(reproject(source, target, dim, [val]))
end
function _reproject(source::GeoFormat, target::GeoFormat, dim, vals::NTuple{N}) where N
    reps = reproject(source, target, dim, [vals...])
    return ntuple(x -> reps[x], N)
end
function _reproject(source::GeoFormat, target::GeoFormat, dim::Union{XDim,YDim}, vals::AbstractArray) 
    reshape(reproject(source, target, dim, vec(vals)), size(vals))
end
function _reproject(source::Nothing, target::GeoFormat, dim::Union{XDim,YDim}, vals::AbstractArray) 
    reshape(reproject(source, target, dim, vec(vals)), size(vals))
end
function _reproject(source::GeoFormat, target::GeoFormat, dim::Union{XDim,YDim}, vals::AbstractVector) 
    error("Rasters.jl requires backends to be loaded externally as of version 0.8. Run `using ArchGDAL` to use `reproject`")
end

# Guess the step for arrays
stepof(A::AbstractArray) = (last(A) - first(A)) / (length(A) - 1)
stepof(A::AbstractRange) = step(A)

isregular(A::AbstractRange) = true
function isregular(A::AbstractArray)
    step = stepof(A)
    for i in eachindex(A)[2:end]
        if !(A[i] - A[i-1] ≈ step)
            return false 
        end
    end
    return true
end
