"""
    points(A::AbstractRaster; dims=(YDim, XDim), ignore_missing) => Array{Tuple}

Returns a generator of the points in `A` for dimensions in `dims`,
where points are a tuple of the values in each specified dimension
index.

# Keywords

- `dims` the dimensions to return points from. The first slice of other
    layers will be used.
- `ignore_missing`: wether to ignore missing values in the array when considering
    points. If `true`, all points in the dimensions will be returned, if `false`
    only the points that are not `=== missingval(A)` will be returned.

The order of `dims` determines the order of the points.

$EXPERIMENTAL
"""
function points(A::AbstractRaster; ignore_missing=false, order=(XDim, YDim, ZDim))
    ignore_missing ? _points(A; order) : _points_missing(A; order)
end
function points(dims::DimTuple; order=(XDim, YDim, ZDim))
    return DimPoints(dims; order=DD.dims(dims, order))
end

_points(A::AbstractRaster; kw...) = points(dims(A); kw...)
function _points_missing(A::AbstractRaster; order)
    points = DimPoints(A; order=dims(A, order))
    function ordered_point_or_missing(I) 
        isequal(A[I], missingval(A)) ? missing : points[I]
    end
    return (ordered_point_or_missing(I) for I in CartesianIndices(A))
end
