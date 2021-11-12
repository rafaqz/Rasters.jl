"""
    reproject(source::GeoFormat, target::GeoFormat, dim::Dimension, val)

`reproject` uses ArchGDAL.reproject, but implemented for a reprojecting
a value array of values, a single dimension at a time.
"""
reproject(source, target, dim, val) = val
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
