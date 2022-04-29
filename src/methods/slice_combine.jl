"""
    slice(A::Union{AbstractRaster,AbstractRasterStack,AbstracRasterSeries}, dims) => RasterSeries

Slice an object along some dimension/s, lazily using `view`.

For a single `Raster` or `RasterStack` this will return a `RasterSeries` of
`Raster` or `RasterStack` that are slices along the specified dimensions.

For a `RasterSeries`, the output is another series where the child objects are sliced and the
series dimensions index is now of the child dimensions combined. `slice` on a `RasterSeries`
with no dimensions will slice along the dimensions shared by both the series and child object.

$EXPERIMENTAL
"""
slice(x::RasterStackOrArray, dims) = slice(x, (dims,))
# Slice an array or stack into a series
function slice(x::RasterStackOrArray, dims::Tuple)
    # Make sure all dimensions in `dims` are in `x`
    all(hasdim(x, dims)) || _errordimsnotfound(dims, DD.dims(x))
    # Define dimensions and data for the sliced RasterSeries
    seriesdims = DD.dims(x, dims)
    # series data is a generator of view slices
    seriesdata = map(DimIndices(seriesdims)) do ds
        view(x, ds...)
    end
    return RasterSeries(seriesdata, seriesdims)
end
# Slice an existing series into smaller slices
slice(ser::AbstractRasterSeries, dims) = cat(map(x -> slice(x, dims), ser)...; dims=dims)

@noinline _errordimsnotfound(targets, dims) =
    throw(ArgumentError("Dimensions $(map(DD.dim2key, targets)) were not found in $(map(DD.dim2key, dims))"))

# By default, combine all the RasterSeries dimensions and return a Raster or RasterStack
combine(ser::AbstractRasterSeries) = combine(ser, dims(ser))
# Fold over all the dimensions, combining the series one dimension at a time
combine(ser::AbstractRasterSeries, dims::Tuple) = foldl(combine, dims; init=ser)
# Slice the N-dimensional series into an array of 1-dimensional
# series, and combine them, returning a new series with 1 less dimension.
function combine(ser::AbstractRasterSeries{<:Any,M}, dim::Union{Dimension,DD.DimType,Val,Symbol}) where M
    od = otherdims(ser, dim)
    slices = map(d -> view(ser, d...), DimIndices(od))
    newchilren = map(s -> combine(s, dim), slices)
    return rebuild(ser; data=newchilren, dims=od)
end
# Actually combine a 1-dimensional series with `cat`
function combine(ser::AbstractRasterSeries{<:Any,1}, dim::Union{Dimension,DD.DimType,Val,Symbol})
    dim = DD.dims(ser, dim)
    D = DD.basetypeof(dim)
    x = foldl(ser) do acc, x
        # May need to reshape to match acc
        cat(acc, _maybereshape(x, acc, dim); dims=D)
    end
    return set(x, D => dims(ser, dim))
end

function _maybereshape(A::AbstractRaster{<:Any,N}, acc, dim) where N
    if ndims(acc) != ndims(A)
        newdata = reshape(parent(A), Val{N+1}())
        d = if hasdim(refdims(A), dim)
            dims(refdims(A), dim)
        else
            DD.basetypeof(dim)(1:1; lookup=NoLookup())
        end
        newdims = (DD.dims(A)..., d)
        return rebuild(A; data=newdata, dims=newdims)
    else
        return A
    end
end
function _maybereshape(st::AbstractRasterStack, acc, dim)
    map((s, a) -> _maybereshape(s, a, dim), st, acc)
end

# See iterate(::GridChunks) in Diskarrays.jl
function _chunk_inds(g, ichunk)
    outinds = map(ichunk.I, g.chunksize, g.parentsize, g.offset) do ic, cs, ps, of
        max((ic - 1) * cs + 1 -of, 1):min(ic * cs - of, ps)
    end
end
