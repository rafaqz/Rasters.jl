
function Rasters._reproject(source::GeoFormat, target::GeoFormat, ::XDim, vals::AbstractVector) 
    paired_vals = map(v -> (v, zero(v)), vals)
    push!(paired_vals, (first(vals), one(first(vals))))
    reps = AG.reproject(paired_vals, source, target; order=:trad)
    at_one = pop!(reps)
    first(reps)[1] == at_one[1] || _reproject_crs_error(source, target)
    return map(r -> r[1], reps)
end
function Rasters._reproject(source::GeoFormat, target::GeoFormat, dim::YDim, vals::AbstractVector)
    paired_vals = map(v -> (zero(v), v), vals)
    push!(paired_vals, (one(first(vals)), first(vals)))
    reps = AG.reproject(paired_vals, source, target; order=:trad)
    at_one = pop!(reps)
    first(reps)[2] == at_one[2] || _reproject_crs_error(source, target)
    return map(r -> r[2], reps)
end

_reproject_crs_error(source, target) = throw(ArgumentError("Cannot reproject from: \n $source \nto: \n $target")) 
