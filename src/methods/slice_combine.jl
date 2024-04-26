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
    combine(A::Union{AbstractRaster,AbstractRasterStack,AbstracRasterSeries}, [dims]) => Raster

Combine a `RasterSeries` along some dimension/s, creating a new `Raster` or `RasterStack`,
depending on the contents of the series.

If `dims` are passed, only the specified dimensions will be combined
with a `RasterSeries` returned, unless `dims` is all the dims in the series.

$EXPERIMENTAL
"""
function combine(ser::AbstractRasterSeries, dims)
    ods = otherdims(ser, dims)
    if length(ods) > 0
        map(DimIndices(ods)) do D
            combine(view(ser, D...))
        end |> RasterSeries
    else
        combine(ser)
    end
end
function combine(ser::AbstractRasterSeries{<:Any,N}) where N
    r1 = first(ser)
    dest = _alloc_combine_dest(ser)
    for sD in DimIndices(ser)
        rD = map(d -> rebuild(d, :), DD.dims(r1))
        source = ser[sD...]
        if dest isa RasterStack
            foreach(layers(source), layers(dest)) do source_r, dest_r
                view(dest_r, rD..., sD...) .= source_r
            end
        else
            # TODO we shouldn't need a view here??
            view(dest, rD..., sD...) .= source
        end
    end
    return dest
end

_alloc_combine_dest(s::AbstractRasterSeries) = _alloc_combine_dest(first(s), s)
_alloc_combine_dest(r::AbstractRaster, s) = similar(r, (dims(r)..., dims(s)...))
_alloc_combine_dest(r::AbstractRasterStack, s) = map(r -> _alloc_combine_dest(r, s), first(s))

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
