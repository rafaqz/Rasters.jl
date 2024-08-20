"""
    slice(A::Union{AbstractRaster,AbstractRasterStack,AbstracRasterSeries}, dims) => RasterSeries

Slice views along some dimension/s to obtain a `RasterSeries` of the slices.

For a `Raster` or `RasterStack` this will return a `RasterSeries` of
`Raster` or `RasterStack` that are slices along the specified dimensions.

For a `RasterSeries`, the output is another series where the child objects are sliced and the
series dimensions index is now of the child dimensions combined. `slice` on a `RasterSeries`
with no dimensions will slice along the dimensions shared by both the series and child object.

$EXPERIMENTAL
"""
slice(x::RasterStackOrArray, dim) = slice(x, (dim,))
function slice(x::RasterStackOrArray, dims::Tuple)
    seriesdims = DD.dims(x, dims)
    seriesdata = eachslice(x; dims)
    @static if VERSION >= v"1.9"
        if seriesdata isa Slices
            return RasterSeries(seriesdata, seriesdims)
        end
    end
    return RasterSeries(collect(seriesdata), seriesdims)
end
# Slice an existing series into smaller slices
slice(ser::AbstractRasterSeries, dims) = cat(map(x -> slice(x, dims), ser)...; dims=dims)

@noinline _errordimsnotfound(targets, dims) =
    throw(ArgumentError("Dimensions $(map(name, targets)) were not found in $(map(name, dims))"))

"""
    combine(A::AbstracRasterSeries; [dims], [lazy]) => Raster

Combine a `RasterSeries` along some dimension/s, creating a new `Raster` or `RasterStack`,
depending on the contents of the series.

If `dims` are passed, only the specified dimensions will be combined
with a `RasterSeries` returned, unless `dims` is all the dims in the series.

If `lazy`, concatenate lazily. The default is to concatenate lazily for lazy `Raster`s and eagerly otherwise.

$EXPERIMENTAL
"""
combine(ser::AbstractRasterSeries, dims; kw...) = combine(ser; dims, kw...)
function combine(ser::AbstractRasterSeries; dims, lazy = isdisk(ser))
    ods = otherdims(ser, dims)
    if length(ods) > 0
        map(x -> combine(x; lazy), eachslice(ser; dims = ods))
    else
        combine(ser; lazy)
    end
end
function combine(ser::AbstractRasterSeries; lazy = isdisk(ser))
    ras1 = first(ser)
    alldims = (dims(ras1)..., dims(ser)...)
    ser_res = DD._insert_length_one_dims(ser, alldims)
    data = DA.ConcatDiskArray(ser_res)
    data = lazy ? data : collect(data)
    return rebuild(ras1; data, dims = alldims)
end
