
const DimOrTuple = Union{AbstractDimension,Tuple{Vararg{<:AbstractDimension}}}
const IntOrTuple = Union{Int,Tuple{Vararg{<:Int}}}

"""
    aggregate(x, method, scale)

Aggregate array, or the arrays in a stack or series by some scale,
using some aggregation function or a [`Locus`] type to specify the position
to sample from.
"""
function aggregate end

aggregate(x; method=Center(), scale) = aggregate(x, method, scale)

aggregate(series::AbstractGeoSeries, method, scale::DimOrTuple) =
    map(x -> aggregate(x, method, scale), series)
aggregate(series::AbstractGeoSeries, method, scale) =
    map(x -> aggregate(x, method, scale), series)
aggregate(stack::AbstractGeoStack, method, scale) = begin
    data = map(NamedTuple{keys(stack)}(keys(stack))) do key
        aggregate(stack[key], method, scale)
    end
    GeoStack(stack; data=data)
end

aggregate(src::AbstractGeoArray, method, scale::DimOrTuple) =
    aggregate(src, method, dims2indices(src, scale))
aggregate(src::AbstractGeoArray, method, scale::IntOrTuple) =
    aggregate(GeoArray(src), method, scale)
aggregate(src::GeoArray, method, scale::IntOrTuple) =
    aggregate!(init_aggregation(src, method, scale), src, method, scale)

"""
    downsample!(out::AbstractMatrix, a::AbstractMatrix, method, scale)

Downsample matrix `a` to another matrix `out` of the correct size.

- `method` is a function such as mean or sum that can combine the
    value of multiple cells to generate the aggregated cell.
- `scale` is the aggregation factor.
"""
aggregate!(dst, src, method, scale::DimOrTuple) =
    aggregate!(dst, src, method, dims2indices(src, scale))
# Loci methods
aggregate!(dst::AbstractGeoArray, src, locus::Locus, scale) =
    aggregate!(dst, src, (locus,), scale)
aggregate!(dst::AbstractGeoArray, src, loci::Tuple{Locus,Vararg}, scale) = begin
    offsets = beginoffset.(loci, scale)
    for I in CartesianIndices(dst)
        dst[I] = src[(upsample.(Tuple(I), scale) .+ offsets)...]
    end
    dst
end
# Func/functor methods
aggregate!(dst::AbstractGeoArray, src, f, scale) = begin
    for I in CartesianIndices(dst)
        topleft = upsample.(Tuple(I), scale)
        bottomright = topleft .+ scale .- 1
        dst[I] = f(view(src, map(:, topleft, bottomright)...))
    end
    dst
end

"""
    initaggregation(A, scale)

Generate an array for aggregating array `A` by `scale`.
"""
init_aggregation(A::AbstractArray, method, scale) =
    init_aggregation(A::AbstractArray, (method,), scale)
init_aggregation(A::AbstractArray, method::Tuple, scale) = begin
    _dims = aggregate.(dims(A), method, scale)
    # Dim aggregation determines the array size 
    lengths = map(length, _dims)
    _data = similar(data(A), lengths)
    rebuild(A; data=_data, dims=_dims)
end

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

aggregate(dim::AbstractDimension, method, scale) =
    aggregate(grid(dim), dim, method, scale)
aggregate(grid, dim::AbstractDimension, method, scale) = begin
    start, stop = endpoints(dim, method, scale)
    rebuild(dim, val(dim)[start:scale:stop])
end

endpoints(dim, method, scale) = begin
    start = firstindex(dim) + beginoffset(dim, method, scale)
    stop = (length(dim) ÷ scale) * scale
    start, stop
end

beginoffset(dim, locus::Locus, scale) = beginoffset(locus, scale)
beginoffset(dim, method, scale) = beginoffset(locus(dim), scale)
beginoffset(locus::Start, scale) = 0
beginoffset(locus::End, scale) = scale - 1
beginoffset(locus::Center, scale) = scale ÷ 2
