function _geotransform2affine(gt::AbstractVector)
    M = [gt[GDAL_WE_RES] gt[GDAL_ROT1]; gt[GDAL_ROT2] gt[GDAL_NS_RES]]
    v = [gt[GDAL_TOPLEFT_X], gt[GDAL_TOPLEFT_Y]]
    return CoordinateTransformations.AffineMap(M, v)
end

function _affine2geotransform(am::CoordinateTransformations.AffineMap)
    M = am.linear
    v = am.translation
    gt = zeros(6)
    gt[GDAL_TOPLEFT_X] = v[1]
    gt[GDAL_WE_RES] = M[1, 1]
    gt[GDAL_ROT1] = M[1, 2]
    gt[GDAL_TOPLEFT_Y] = v[2]
    gt[GDAL_ROT2] = M[2, 1]
    gt[GDAL_NS_RES] = M[2, 2]
    return gt
end

function _affine_extrema(am::CoordinateTransformations.AffineMap, lookup_x, lookup_y)
    # Not 100% sure this holds in all cases
    minx = first(lookup_x) - 1
    maxx = last(lookup_x)
    miny = first(lookup_y) - 1
    maxy = last(lookup_y)
    extrema = am((minx, miny)), am((maxx, maxy)), am((minx, maxy)), am((maxx, miny))
end
