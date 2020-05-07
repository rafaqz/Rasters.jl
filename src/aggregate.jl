
const DimOrDimTuple = Union{Dimension,Tuple{Vararg{<:Dimension}}}
const IntOrIntTuple = Union{Int,Tuple{Vararg{<:Int}}}

"""
    aggregate(x, method, scale)

Aggregate array, or all arrays in a stack or series, by some scale.
This uses a `Array` aggregation function like `mean`, or a [`Locus`] type to
specify a single position to sample from. Return values are `GeoArray`,
`GeoStack` or `GeoSeries` depending on the type of `x`.

- `method` is a function such as mean or sum that can combine the
    value of multiple cells to generate the aggregated cell, or a loci
    like `Start` or `Center()` that species where to sample from in the interval.
- `scale` is the aggregation factor, which can be an integer, a tuple of integers
  for each dimension, or any `Dimension`, `Selector` or `Int` combination you can
  usually use in `getindex`. Using a `Selector` will determine the scale by the
  distance from the start of the index.
"""
function aggregate end

"""
    aggregate(method, series::AbstractGeoSeries, scale)

Aggregate an AbstractGeoSeries
"""
aggregate(method, series::AbstractGeoSeries, scale, args...) =
    rebuild(series, [aggregate(method, series[i], scale, args...) for i in 1:length(series)])
"""
    aggregate(method, stack::AbstractGeoStack, scale)

Aggregate an AbstractGeoStack
"""
aggregate(method, stack::AbstractGeoStack, scale, keys=keys(stack)) = begin
    data = map(NamedTuple{keys}(keys)) do key
        aggregate(method, stack[key], scale)
    end
    GeoStack(stack; data=data)
end
aggregate(method, src::DiskGeoArray, scale) =
    aggregate(method, GeoArray(src), scale)

# DimensionalData methods
"""
    aggregate(method, src::AbstractDimensionalArray, scale)

Aggregate an AbstractDimensionalArray
"""
aggregate(method, src::AbstractDimensionalArray, scale) =
    aggregate!(method, ag_array(method, src, scale), src, scale)
"""
    aggregate(method, dim::Dimension, scale)

Aggregate a Dimension
"""
aggregate(method, dim::Dimension, scale) = begin
    start, stop = endpoints(dim, method, scale)
    rebuild(dim, val(dim)[start:scale:stop], aggregate(method, mode(dim), scale))
end

endpoints(dim, method, scale) = begin
    start = firstindex(dim) + beginoffset(method, dim, scale)
    stop = (length(dim) ÷ scale) * scale
    start, stop
end

"""
    aggregate(method, dim::IndexMode, scale)

Aggregate an IndexMode
"""
aggregate(method, mode::IndexMode, scale) = mode 
aggregate(method, mode::AbstractSampled, scale) =
    rebuild(mode; span=aggregate(method, span(mode), scale)) 

"""
    aggregate(method, dim::Span, scale)

Aggregate a Span
"""
aggregate(method, span::Span, scale) = span
aggregate(method, span::Regular, scale) = Regular(val(span) * scale) 

"""
    aggregate!(dst::AbstractDimensionalArray, src::AbstractDimensionalArray, method, scale)

Aggregate array `src` to array `dst` by some scale.
This uses an aggregation function like `mean` or a [`Locus`] type to
specify a position to sample from.

- `method` is a function such as mean or sum that can combine the
    value of multiple cells to generate the aggregated cell, or a loci
    like `Start` or `Center()` that species where to sample from in the interval.
- `scale` is the aggregation factor, which can be an integer, or a tuple of an
  `Dimension`, `Selector` or `Int` combination you can usually use in `getindex`.
  Using a `Selector` will determine the scale by the distance from the start of the index.
"""
aggregate!(locus::Locus, dst::AbstractDimensionalArray, src, scale) =
    aggregate!((locus,), dst, src, scale)
aggregate!(loci::Tuple{Locus,Vararg}, dst::AbstractDimensionalArray, src, scale) = begin
    intscale = scale2int(dims(src), scale)
    offsets = beginoffset.(loci, intscale)
    for I in CartesianIndices(dst)
        dst[I] = src[(upsample.(Tuple(I), intscale) .+ offsets)...]
    end
    dst
end
# Function/functor methods
aggregate!(f, dst::AbstractDimensionalArray, src, scale) = begin
    intscale = scale2int(dims(src), scale)
    for I in CartesianIndices(dst)
        topleft = upsample.(Tuple(I), intscale)
        bottomright = topleft .+ intscale .- 1
        dst[I] = f(view(src, map(:, topleft, bottomright)...))
    end
    dst
end

# Allocate an array of the correct size to aggregate `A` by `scale`
ag_array(method, A::AbstractDimensionalArray, scale) =
    ag_array((method,), A, scale)
ag_array(method::Tuple, A::AbstractDimensionalArray, scale) = begin
    intscale = scale2int(dims(A), scale)
    # Aggregate the dimensions
    dims_ = aggregate.(method, dims(A), intscale)
    # Dim aggregation determines the array size
    data_ = similar(data(A), map(length, dims_)...)
    rebuild(A; data=data_, dims=dims_)
end

# Convert scale or tuple of scale to integer using dims2indices
scale2int(dims::Tuple, scale::Tuple) = dims2indices(dims, scale)
scale2int(dims::Tuple, scale::Int) = scale

"""
    upsample(index::Int, scale::Int)

Convert indicies from the aggregated array to the larger original array.
"""
upsample(index::Int, scale::Int) = (index - 1) * scale + 1

"""
    downsample(index::Int, scale::Int)

Convert indicies from the original array to the aggregated array.
"""
downsample(index::Int, scale::Int) = (index - 1) ÷ scale + 1

beginoffset(locus::Locus, dim::Dimension, scale) = beginoffset(locus, scale)
beginoffset(method, dim::Dimension, scale) = beginoffset(locus(dim), scale)
beginoffset(locus::Start, scale) = 0
beginoffset(locus::End, scale) = scale - 1
beginoffset(locus::Center, scale) = scale ÷ 2
