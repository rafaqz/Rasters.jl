
Rasters._reproject(source::GeoFormat, target::GeoFormat, dim::XDim, vals::AbstractVector) = _reproject!(source, target, dim, Float64.(vals))
Rasters._reproject(source::GeoFormat, target::GeoFormat, dim::YDim, vals::AbstractVector) = _reproject!(source, target, dim, Float64.(vals))

function _reproject!(source::GeoFormat, target::GeoFormat, ::XDim, vals::Vector{Float64}) 
    # First, construct a transformation from `source` to `target`.
    # TODO: add area of use, so that the transformation is more accurate.
    trans = Proj.Transformation(source, target; always_xy = true, direction = Proj.PJ_FWD)
    # Check that the reprojection is valid, i.e., trans((x, 0))[1] == trans((x, 1))[1]
    first_x = first(vals)
    zero_y = (first_x, zero(first_x))
    one_y = (first_x, one(first_x))
    # Check that the reprojection is valid, i.e., trans((x, 0))[1] == trans((x, 1))[1]
    trans(zero_y)[1] == trans(one_y)[1] || _reproject_crs_error(source, target)
    # Now that we've proved that the reprojection is valid, we can proceed.
    # Here, for efficiency, and to avoid allocations, we'll use `proj_trans_generic`,
    # which mutates values in-place.
    Proj.proj_trans_generic(
        trans.pj,                                    # pointer to transformation
        Proj.PJ_FWD,                                 # direction of transformation
        vals, sizeof(typeof(first_x)), length(vals), # input x values
        C_NULL, 0, 0,                                # default (assumed 0) y values
        C_NULL, 0, 0,                                # default (assumed 0) z values
        C_NULL, 0, 0                                 # default (assumed 0) t values
    )
    return vals
end
function _reproject!(source::GeoFormat, target::GeoFormat, dim::YDim, vals::AbstractVector)
    # First, construct a transformation from `source` to `target`.
    # TODO: add area of use, so that the transformation is more accurate.
    trans = Proj.Transformation(source, target; always_xy = true, direction = Proj.PJ_FWD)
    # Check that the reprojection is valid, i.e., trans((0, y))[2] == trans((1, y))[2]
    first_y = first(vals)
    zero_x = (zero(first_y), first_y)
    one_x = (one(first_y), first_y)
    # Check that the reprojection is valid, i.e., trans((0, y))[2] == trans((1, y))[2]
    trans(zero_x)[2] == trans(one_x)[2] || _reproject_crs_error(source, target)
    # Now that we've proved that the reprojection is valid, we can proceed.
    # Here, for efficiency, and to avoid allocations, we'll use `proj_trans_generic`,
    # which mutates values in-place.
    Proj.proj_trans_generic(
        trans.pj,                                    # pointer to transformation
        Proj.PJ_FWD,                                 # direction of transformation
        C_NULL, 0, 0,                                # default (assumed 0) x values
        vals, sizeof(typeof(first_y)), length(vals), # input y values
        C_NULL, 0, 0,                                # default (assumed 0) z values
        C_NULL, 0, 0                                 # default (assumed 0) t values
    )
    return vals
end

_reproject_crs_error(source, target) = 
    throw(ArgumentError("Cannot reproject from: \n $source \nto: \n $target, coordinate reference systems are not aligned on all axes. You may need to use `resample` instead")) 
