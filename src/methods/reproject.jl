"""
    reproject(target::GeoFormat, x)

Reproject the dimensions of `x` to a different crs.

# Arguments

- `target`: any crs in a GeoFormatTypes.jl wrapper, e.g. `EPSG`, `WellKnownText`, `ProjString`.
- `x`: a `Dimension`, `Tuple` of `Dimension`, `Raster` or `RasterStack`.

Dimensions without an `AbstractProjected` lookup (such as a `Ti` dimension)
are silently returned without modification.
"""
reproject(target::GeoFormat, x) = rebuild(x; dims=reproject(target, dims(x)))
reproject(target::GeoFormat, dims::Tuple) = map(d -> reproject(target, d), dims)
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
reproject(target::GeoFormat, l::LookupArray) = l

"""
    reproject(source::GeoFormat, target::GeoFormat, dim::Dimension, val)

`reproject` uses ArchGDAL.reproject, but implemented for a reprojecting
a value array of values, a single dimension at a time.
"""
reproject(source::GeoFormat, target::GeoFormat, dim, val) = val
reproject(source::Nothing, target::GeoFormat, dim, val) = val
reproject(source::GeoFormat, target::Nothing, dim, val) = val
reproject(source::Nothing, target::Nothing, dim, val) = val
function reproject(source::GeoFormat, target::GeoFormat, dim::Union{XDim,YDim}, val::Number)
    # This is a dumb way to do this. But it save having to inspect crs, 
    # and prevents reprojections that dont make sense from working.
    # A better method for this should be implemented in future.
    return first(reproject(source, target, dim, [val]))
end
function reproject(source::GeoFormat, target::GeoFormat, dim, vals::NTuple{N}) where N
    reps = reproject(source, target, dim, [vals...])
    return ntuple(x -> reps[x], N)
end
function reproject(source::GeoFormat, target::GeoFormat, dim::Union{XDim,YDim}, vals::AbstractArray) 
    reshape(reproject(source, target, dim, vec(vals)), size(vals))
end
function reproject(source::Nothing, target::GeoFormat, dim::Union{XDim,YDim}, vals::AbstractArray) 
    reshape(reproject(source, target, dim, vec(vals)), size(vals))
end
function reproject(source::GeoFormat, target::GeoFormat, ::XDim, vals::AbstractVector) 
    paired_vals = map(v -> (v, zero(v)), vals)
    push!(paired_vals, (first(vals), one(first(vals))))
    reps = AG.reproject(paired_vals, source, target; order=:trad)
    at_one = pop!(reps)
    first(reps)[1] == at_one[1] || _rep_error(source, target)
    return map(r -> r[1], reps)
end
function reproject(source::GeoFormat, target::GeoFormat, dim::YDim, vals::AbstractVector)
    paired_vals = map(v -> (zero(v), v), vals)
    push!(paired_vals, (one(first(vals)), first(vals)))
    reps = AG.reproject(paired_vals, source, target; order=:trad)
    at_one = pop!(reps)
    first(reps)[2] == at_one[2] || _rep_error(source, target)
    return map(r -> r[2], reps)
end

_rep_error(source, target) = throw(ArgumentError("Cannot reproject from: \n $source \nto: \n $target")) 

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
