
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
aggregate(method, series::AbstractGeoSeries, scale, args...; progress=true, kwargs...) = begin
    f = i -> aggregate(method, series[i], scale, args...; progress=false, kwargs...)
    data = if progress
        @showprogress "Agregating series..." map(f, 1:length(series))
    else
        map(f, 1:length(series))
    end
    rebuild(series, data)
end
"""
    aggregate(method, stack::AbstractGeoStack, scale)

Aggregate an AbstractGeoStack
"""
aggregate(method, stack::AbstractGeoStack, scale;
          keys=keys(stack), progress=true) = begin
    f = key -> aggregate(method, stack[key], scale)
    keys_nt = NamedTuple{keys}(keys)
    data = if progress
        @showprogress "Agregating stack..." map(f, keys_nt)
    else
        map(f, keys_nt)
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
    aggregate!(method, alloc_ag(method, src, scale), src, scale)
"""
    aggregate(method, dim::Dimension, scale)

Aggregate a Dimension
"""
aggregate(method, dim::Dimension, scale) = begin
    intscale = scale2int(dim, scale)
    start, stop = endpoints(dim, method, intscale)
    rebuild(dim, val(dim)[start:scale:stop], aggregate(method, mode(dim), intscale))
end

endpoints(dim, method, scale) = begin
    start = firstindex(dim) + agoffset(method, dim, scale)
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
    offsets = agoffset.(loci, intscale)
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
alloc_ag(method, A::AbstractDimensionalArray, scale) =
    alloc_ag((method,), A, scale)
alloc_ag(method::Tuple, A::AbstractDimensionalArray, scale) = begin
    intscale = scale2int(dims(A), scale)
    # Aggregate the dimensions
    dims_ = aggregate.(method, dims(A), intscale)
    # Dim aggregation determines the array size
    data_ = similar(data(A), map(length, dims_)...)
    rebuild(A; data=data_, dims=dims_)
end

# Convert scale or tuple of scale to integer using dims2indices
scale2int(dims, scale::Tuple) = dims2indices(dims, scale)
scale2int(dims, scale::Int) = scale


"""
    disaggregate(x, method, scale)

Disagregate array, or all arrays in a stack or series, by some scale.
This uses a `Array` aggregation function like `mean`, or a [`Locus`] type to
specify a single position to sample from. Return values are `GeoArray`,
`GeoStack` or `GeoSeries` depending on the type of `x`.

- `method` is a function such as mean or sum that can combine the
    value of multiple cells to generate the disaggregated cell, or a loci
    like `Start` or `Center()` that species where to sample from in the interval.
- `scale` is the aggregation factor, which can be an integer, a tuple of integers
  for each dimension, or any `Dimension`, `Selector` or `Int` combination you can
  usually use in `getindex`. Using a `Selector` will determine the scale by the
  distance from the start of the index.
"""
function disaggregate end

"""
    disaggregate(method, series::AbstractGeoSeries, scale)

Disagregate an AbstractGeoSeries
"""
disaggregate(method, series::AbstractGeoSeries, scale; progress=true, kwargs...) = begin
    f = i -> disaggregate(method, series[i], scale; progress=false, kwargs...)
    data = if progress
        @showprogress "Disgregating series..." map(f, 1:length(series))
    else
        map(f, 1:length(series))
    end
end
"""
    disaggregate(method, stack::AbstractGeoStack, scale)

Disagregate an AbstractGeoStack
"""
disaggregate(method, stack::AbstractGeoStack, scale; keys=keys(stack), progress=true) = begin
    f = key -> disaggregate(method, stack[key], scale)
    keys_nt = NamedTuple{keys}(keys)
    data = if progress
        @showprogress "Disagregating stack..." map(f, keys_nt)
    else
        map(f, keys_nt)
    end
    GeoStack(stack; data=data)
end
disaggregate(method, src::DiskGeoArray, scale) =
    disaggregate(method, GeoArray(src), scale)

# DimensionalData methods
"""
    disaggregate(method, src::AbstractDimensionalArray, scale)

Disagregate an AbstractDimensionalArray
"""
disaggregate(method, src::AbstractDimensionalArray, scale) =
    disaggregate!(method, alloc_disag(method, src, scale), src, scale)

"""
    disaggregate(method, dim::Dimension, scale)

Disagregate a Dimension
"""
disaggregate(locus::Locus, dim::Dimension, scale) = begin
    intscale = scale2int(dim, scale)
    len = length(dim) * intscale
    step_ = step(mode(dim)) / intscale
    start = dim[1] - agoffset(locus, intscale) * step_
    stop = start + (len - 1)  * step_
    index = LinRange(start, stop, len)
    rebuild(dim, index, disaggregate(locus, mode(dim), intscale))
end
disag_index(locus::Start, dim, scale) =
    LinRange(dim[1], dim[end] + (scale - 1) * step(dim), length(dim) * scale)
disag_index(locus::End, dim, scale) =
    LinRange(dim[1], dim[end] + (scale - agoffset(locus)) * step(dim), length(dim) * scale)



"""
    disaggregate(method, dim::IndexMode, scale)

Disagregate an IndexMode
"""
disaggregate(method, mode::IndexMode, scale) = mode
disaggregate(method, mode::AbstractSampled, scale) =
    rebuild(mode; span=disaggregate(method, span(mode), scale))

"""
    disaggregate(method, dim::Span, scale)

Disagregate a Span
"""
disaggregate(method, span::Span, scale) = span
disaggregate(method, span::Regular, scale) = Regular(val(span) / scale)

"""
    disaggregate!(dst::AbstractDimensionalArray, src::AbstractDimensionalArray, method, scale)

Disagregate array `src` to array `dst` by some scale.
This uses an aggregation function like `mean` or a [`Locus`] type to
specify a position to sample from.

- `method` is a function such as mean or sum that can combine the
    value of multiple cells to generate the disaggregated cell, or a loci
    like `Start` or `Center()` that species where to sample from in the interval.
- `scale` is the aggregation factor, which can be an integer, or a tuple of an
  `Dimension`, `Selector` or `Int` combination you can usually use in `getindex`.
  Using a `Selector` will determine the scale by the distance from the start of the index.
"""
disaggregate!(locus::Locus, dst::AbstractDimensionalArray, src, scale) =
    disaggregate!((locus,), dst, src, scale)
disaggregate!(loci::Tuple{Locus,Vararg}, dst::AbstractDimensionalArray, src, scale) = begin
    intscale = scale2int(dims(src), scale)
    for I in CartesianIndices(dst)
        dst[I] = src[(downsample.(Tuple(I), intscale))...]
    end
    dst
end

# Allocate an array of the correct size to aggregate `A` by `scale`
alloc_disag(method, A::AbstractDimensionalArray, scale) =
    alloc_disag((method,), A, scale)
alloc_disag(method::Tuple, A::AbstractDimensionalArray, scale) = begin
    intscale = scale2int(dims(A), scale)
    dims_ = disaggregate.(method, dims(A), intscale)
    # Dim aggregation determines the array size
    data_ = similar(data(A), map(length, dims_)...)
    rebuild(A; data=data_, dims=dims_)
end


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

agoffset(locus::Locus, dim::Dimension, scale) = agoffset(locus, scale)
agoffset(method, dim::Dimension, scale) = agoffset(locus(dim), scale)
agoffset(locus::Start, scale) = 0
agoffset(locus::End, scale) = scale - 1
agoffset(locus::Center, scale) = scale ÷ 2
