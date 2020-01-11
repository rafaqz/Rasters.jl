
"""
    aggregate(x, aggregator, scale)

Aggregate array, or the arrays in a stack or series by some scale,
using some aggregation function or a [`Locus`] type to specify the position
to sample from.
"""
function aggregate end

aggregate(stack::AbstractGeoStack, aggregator, scale) = begin
    data = map(NamedTuple{keys(stack)}(keys(stack))) do key
        aggregate(GeoArray(stack[key]), aggregator, scale)
    end
    GeoStack(stack; data=data)
end

const DimOrTuple = Union{AbstractDimension,Tuple{Vararg{<:AbstractDimension}}}

aggregate(src::AbstractGeoArray, aggregator, scale::DimOrTuple) =
    aggregate(src, aggregator, dims2indices(src, scale))
aggregate(src::AbstractGeoArray, aggregator, scale) =
    aggregate!(init_aggregation(src, aggregator, scale), src, aggregator, scale)

"""
    downsample!(out::AbstractMatrix, a::AbstractMatrix, aggregator, scale)

Downsample matrix `a` to another matrix `out` of the correct size.

- `aggregator` is a function such as mean or sum that can combine the
    value of multiple cells to generate the aggregated cell.
- `scale` is the aggregation factor.
"""
aggregate!(dst, src, aggregator, scale::DimOrTuple) =
    aggregate!(dst, src, aggregator, dims2indices(src, scale))
aggregate!(dst::AbstractGeoArray, src, aggregator::Locus, scale) =
    aggregate!(dst::AbstractGeoArray, src, (aggregator,), scale)
aggregate!(dst::AbstractGeoArray, src, aggregator::Tuple{Locus,Vararg}, scale) = begin
    offset = locusoffset.(aggregator, scale)
    for I in CartesianIndices(dst)
        dst[I] = src[(upsample.(Tuple(I), scale) .+ offset)...]
    end
    dst
end

aggregate!(dst::AbstractGeoArray, src, aggregator, scale) = begin
    for I in CartesianIndices(dst)
        I1 = upsample.(Tuple(I), scale)
        I2 = min.(size(src), I1 .+ scale .- 1)
        dst[I] = aggregator(view(src, map(:, I1, I2)...))
    end
    dst
end

"""
    initaggregation(A, scale)

Generate an array for aggregating array `A` by `scale`.
"""
init_aggregation(A::AbstractArray, aggregator, scale) = begin
    _dims = map(aggregate, dims(A), scale)
    rebuild(A; data=similar(A, size(A) .÷ scale), dims=_dims)
end

aggregate(dim::AbstractDimension, scale) =
    rebuild(dim; val=dim[firstindex(dim):scale:lastindex(dim)])

"""
    upsample(index, scale)

Convert indicies from the aggregated array to the larger original array.
"""
upsample(index, scale) = (index - 1) * scale + 1

"""
    downsample(index, scale)

Convert indicies from the original array to the aggregated array.
"""
downsample(index, scale) = (index - 1) ÷ scale + 1


aggregate(dim::AbstractDimension, aggregator, scale) =
    aggregate(grid(dim), dim, aggregator, scale)
aggregate(grid, dim::AbstractDimension, aggregator, scale) =
    aggregate(val(dim), grid, dim, aggregator, scale)

aggregate(index::AbstractRange, grid, dim, aggregator, scale) = begin
    offset = locusoffset(aggregator, scale)
    rebuild(dim, LinRange(index[1 + offset], index[end - scale + 1 + offset], length(index) ÷ scale))
end

aggregate(index::AbstractArray, grid, dim, aggregator, scale) = begin
    offset = locusoffset(aggregator, scale)
    rebuild(dim, [index[(i - 1) * scale + offset] for i in 1:(length(index) ÷ scale)])
end

locusoffset(locus::Start, scale) = 0
locusoffset(locus::End, scale) = scale - 1
locusoffset(locus::Center, scale) = scale ÷ 2
