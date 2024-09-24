
const DimOrDimTuple = Union{Dimension,Tuple{Vararg{<:Dimension}}}
const IntOrIntTuple = Union{Int,Tuple{Vararg{<:Int}}}

struct Ag end
struct DisAg end

"""
    aggregate(method, object, scale; filename, progress, skipmissing)

Aggregate a `Raster`, or all arrays in a `RasterStack` or `RasterSeries`, by `scale` using
`method`.

# Arguments

- `method`: a function such as `mean` or `sum` that can combine the
  value of multiple cells to generate the aggregated cell, or a [`Locus`]($DDlocusdocs)
  like `Start()` or `Center()` that specifies where to sample from in the interval.
- `object`: Object to aggregate, like `AbstractRasterSeries`, `AbstractStack`,
  `AbstractRaster` or `Dimension`.
- `scale`: the aggregation factor, which can be an integer, a tuple of integers
  for each dimension, or any `Dimension`, `Selector` or `Int` combination you can
  usually use in `getindex`. Using a `Selector` will determine the scale by the
  distance from the start of the index.

When the aggregation `scale` of is larger than the array axis, the length of the axis is used.

# Keywords

- `skipmissingval`: if `true`, any `missingval` will be skipped during aggregation, so that
    only areas of all missing values will be aggregated to `missingval`. If `false`, any
    aggregated area containing a `missingval` will be assigned `missingval`.
$FILENAME_KEYWORD
$SUFFIX_KEYWORD
$PROGRESS_KEYWORD

# Example

```jldoctest
using Rasters, RasterDataSources, Statistics, Plots
using Rasters: Center
st = read(RasterStack(WorldClim{Climate}; month=1))
ag = aggregate(Center(), st, (Y(20), X(20)); skipmissingval=true, progress=false)
plot(ag)
savefig("build/aggregate_example.png"); nothing
# output

```

![aggregate](aggregate_example.png)

Note: currently it is faster to aggregate over memory-backed arrays.
Use [`read`](@ref) on `src` before use where required.
"""
function aggregate end
function aggregate(method, series::AbstractRasterSeries, scale, args...;
    progress=true, kw...
)
    f(A) = aggregate(method, A, scale, args...; progress=false, kw...)
    if progress
        ProgressMeter.@showprogress "Aggregating series..." map(f, series)
    else
        map(f, series)
    end
end
function aggregate(method, stack::AbstractRasterStack, scale;
    keys=keys(stack), filename=nothing, suffix=keys, progress=true, kw...
)
    f(A, suffix) = aggregate(method, A, scale; filename, suffix, kw...)

    layers = if progress
        ProgressMeter.@showprogress "Aggregating stack..." map(f, values(stack), Tuple(suffix))
    else
        map(f, values(stack), Tuple(suffix))
    end

    return DD.rebuild_from_arrays(stack, layers)
end
function aggregate(method, src::AbstractRaster, scale;
    suffix=nothing, filename=nothing, kw...
)
    dst = alloc_ag(method, src, scale; filename, suffix, kw...)
    aggregate!(method, dst, src, scale; kw...)
end
aggregate(method, d::Dimension, scale) = rebuild(d, aggregate(method, lookup(d), scale))
function aggregate(method, lookup::Lookup, scale)
    intscale = _scale2int(Ag(), lookup, scale)
    intscale == 1 && return lookup
    start, stop = _endpoints(method, lookup, intscale)
    newlookup = lookup[start:scale:stop]
    if lookup isa AbstractSampled
        sp = aggregate(method, span(lookup), scale)
        return rebuild(newlookup; span=sp)
    else
        return newlookup
    end
end
aggregate(method, span::Span, scale) = span
aggregate(method, span::Regular, scale) = Regular(val(span) * scale)


"""
    aggregate!(method, dst::AbstractRaster, src::AbstractRaster, scale; skipmissingval=false)

Aggregate array `src` to array `dst` by `scale`, using `method`.

# Arguments

- `method`: a function such as `mean` or `sum` that can combine the
  value of multiple cells to generate the aggregated cell, or a [`Locus`]($DDlocusdocs)
  like `Start()` or `Center()` that species where to sample from in the interval.
- `scale`: the aggregation factor, which can be an integer, a tuple of integers
  for each dimension, or any `Dimension`, `Selector` or `Int` combination you can
  usually use in `getindex`. Using a `Selector` will determine the scale by the
  distance from the start of the index in the `src` array.

When the aggregation `scale` of is larger than the array axis, the length of the axis is used.

# Keywords

- `progress`: show a progress bar.
- `skipmissingval`: if `true`, any `missingval` will be skipped during aggregation, so that
    only areas of all missing values will be aggregated to `missingval`. If `false`, any
    aggregated area containing a `missingval` will be assigned `missingval`.

Note: currently it is _much_ faster to aggregate over memory-backed arrays.
Use [`read`](@ref) on `src` before use where required.
"""
function aggregate!(locus::Locus, dst::AbstractRaster, src, scale; kw...)
    aggregate!((locus,), dst, src, scale)
end
function aggregate!(loci::Tuple{Locus,Vararg}, dst::AbstractRaster, src, scale; kw...)
    intscale = _scale2int(Ag(), dims(src), scale)
    offsets = _agoffset.(loci, intscale)
    broadcast!(dst, CartesianIndices(dst)) do I
        val = src[(upsample.(Tuple(I), intscale) .+ offsets)...]
        val === missingval(src) ? missingval(dst) : val
    end
end
# Function/functor methods
function aggregate!(f, dst::AbstractRaster, src, scale; skipmissingval=false)
    intscale = _scale2int(Ag(), dims(src), scale)
    broadcast!(dst, CartesianIndices(dst)) do I
        upper = upsample.(Tuple(I), intscale)
        lower = upper .+ intscale .- 1
        block = if isdisk(src)
            src[map(:, upper, lower)...]
        else
            view(src, map(:, upper, lower)...)
        end
        if skipmissingval
            # All missing values return a missing value
            if all(x -> x === missingval(src), block)
                _missingval_or_missing(dst)
            else
                # Skip missing values
                f((x for x in block if x !== missingval(src)))
            end
        else
            # Any missing values return a missing value
            if any(x -> x === missingval(src), block)
                _missingval_or_missing(dst)
            else
                f(block)
            end
        end
    end
end


"""
    disaggregate(method, object, scale; filename, progress, keys)

Disaggregate array, or all arrays in a stack or series, by some scale.

# Arguments

- `method`: a function such as `mean` or `sum` that can combine the
  value of multiple cells to generate the aggregated cell, or a [`Locus`]($DDlocusdocs)
  like `Start()` or `Center()` that species where to sample from in the interval.
- `object`: Object to aggregate, like `AbstractRasterSeries`, `AbstractStack`,
  `AbstractRaster` or a `Dimension`.
- `scale`: the aggregation factor, which can be an integer, a tuple of integers
  for each dimension, or any `Dimension`, `Selector` or `Int` combination you can
  usually use in `getindex`. Using a `Selector` will determine the scale by the
  distance from the start of the index.

# Keywords

- `progress`: show a progress bar.

Note: currently it is faster to aggregate over memory-backed arrays.
Use [`read`](@ref) on `src` before use where required.

"""
function disaggregate end
function disaggregate(method, series::AbstractRasterSeries, scale; progress=true, kw...)
    f = i -> disaggregate(method, series[i], scale; progress=false, kw...)
    return if progress
        ProgressMeter.@showprogress "Disaggregating series..." map(f, eachindex(series))
    else
        map(f, eachindex(series))
    end
end
function disaggregate(method, stack::AbstractRasterStack, scale;
    keys=keys(stack), suffix=keys, filename=nothing, progress=true
)
    f(A, suffix) = disaggregate(method, A, scale; filename, suffix)

    layers = if progress
        ProgressMeter.@showprogress "Disaggregating stack..." map(f, values(stack), Tuple(suffix))
    else
        map(f, values(stack), Tuple(suffix))
    end
    return DD.rebuild_from_arrays(stack, layers)
end
function disaggregate(method, src::AbstractRaster, scale;
    suffix=nothing, filename=nothing, kw...
)
    dst = alloc_disag(method, src, scale; filename, suffix, kw...)
    disaggregate!(method, dst, src, scale)
end
function disaggregate(locus::Locus, dim::Dimension, scale)
    rebuild(dim, disaggregate(locus, lookup(dim), scale))
end
function disaggregate(locus, lookup::Lookup, scale)
    intscale = _scale2int(DisAg(), lookup, scale)
    intscale == 1 && return lookup

    len = length(lookup) * intscale
    step_ = step(lookup) / intscale
    start = lookup[1] - _agoffset(locus, intscale) * step_
    stop = start + (len - 1)  * step_
    index = LinRange(start, stop, len)
    if lookup isa AbstractSampled
        sp = disaggregate(locus, span(lookup), intscale)
        return rebuild(lookup; data=index, span=sp)
    else
        return rebuild(lookup; data=index)
    end
end

disaggregate(method, span::Span, scale) = span
disaggregate(method, span::Regular, scale) = Regular(val(span) / scale)

"""
    disaggregate!(method, dst::AbstractRaster, src::AbstractRaster, filename, scale)

Disaggregate array `src` to array `dst` by some scale, using `method`.

- `method`: a function such as `mean` or `sum` that can combine the
  value of multiple cells to generate the aggregated cell, or a [`Locus`]($DDlocusdocs)
  like `Start()` or `Center()` that species where to sample from in the interval.
- `scale`: the aggregation factor, which can be an integer, a tuple of integers
  for each dimension, or any `Dimension`, `Selector` or `Int` combination you can
  usually use in `getindex`. Using a `Selector` will determine the scale by the
  distance from the start of the index in the `src` array.

Note: currently it is faster to aggregate over memory-backed arrays.
Use [`read`](@ref) on `src` before use where required.
"""
function disaggregate!(locus::Locus, dst::AbstractRaster, src, scale)
    disaggregate!((locus,), dst, src, scale)
end
function disaggregate!(loci::Tuple{Locus,Vararg}, dst::AbstractRaster, src, scale)
    intscale = _scale2int(DisAg(), dims(src), scale)
    broadcast!(dst, CartesianIndices(dst)) do I
        val = src[(downsample.(Tuple(I), intscale))...]
        val === missingval(src) ? missingval(dst) : val
    end
end

# Allocate an array of the correct size to aggregate `A` by `scale`
alloc_ag(method, A::AbstractRaster, scale; kw...) = alloc_ag((method,), A, scale; kw...)
function alloc_ag(method::Tuple, A::AbstractRaster, scale;
    filename=nothing, suffix=nothing, skipmissingval=nothing
)
    intscale = _scale2int(Ag(), dims(A), scale)
    # Aggregate the dimensions
    dims_ = aggregate.(method, dims(A), intscale)
    # Dim aggregation determines the array size
    sze = map(length, dims_)
    agT = ag_eltype(method, A)
    if missingval(A) isa Nothing
        T = agT
        mv = nothing
    else
        T = promote_type(agT, typeof(missingval(A)))
        mv = convert(T, missingval(A))
    end
    return create(filename, T, dims_; name=name(A), suffix, missingval=mv)
end

# Allocate an array of the correct size to disaggregate `A` by `scale`
function alloc_disag(method, A::AbstractRaster, scale; kw...)
    alloc_disag((method,), A, scale; kw...)
end
function alloc_disag(method::Tuple, A::AbstractRaster, scale;
    filename=nothing, suffix=nothing
)
    intscale = _scale2int(DisAg(), dims(A), scale)
    dims_ = disaggregate.(method, dims(A), intscale)
    # Dim aggregation determines the array size
    sze = map(length, dims_)
    T = ag_eltype(method, A)
    mv = missingval(A) isa Nothing ? nothing : convert(T, missingval(A))
    return create(filename, T, dims_; name=name(A), suffix, missingval=mv)
end

# Handle how methods like `mean` can change the type
ag_eltype(method::Tuple{<:Locus,Vararg}, A) = eltype(A)
function ag_eltype(method::Tuple{<:Any}, A)
    method_returntype = typeof(method[1](zero(eltype(A))))
    promote_type(eltype(A), method_returntype)
end

# Convert indices from the aggregated array to the larger original array.
upsample(index::Int, scale::Int) = (index - 1) * scale + 1
upsample(index::Int, scale::Colon) = index

# Convert indices from the original array to the aggregated array.
downsample(index::Int, scale::Int) = (index - 1) ÷ scale + 1
downsample(index::Int, scale::Colon) = index

# Convert scale or tuple of scale to integer using dims2indices
function _scale2int(x, dims::DimTuple, scale::Tuple)
    map((d, s) -> _scale2int(x, d, s), dims, DD.dims2indices(dims, scale))
end
_scale2int(x, dims::DimTuple, scale::Int) = map(d -> _scale2int(x, d, scale), dims)
_scale2int(x, dims::DimTuple, scale::Colon) = map(d -> _scale2int(x, d, scale), dims)
_scale2int(x, dim::Dimension, scale::Int) = _scale2int(x, lookup(dim), scale)
_scale2int(::Ag, l::Lookup, scale::Int) = scale > length(l) ? length(l) : scale
_scale2int(::DisAg, l::Lookup, scale::Int) = scale
_scale2int(x, dim::Dimension, scale::Colon) = 1

_agoffset(locus::Locus, l::Lookup, scale) = _agoffset(locus, scale)
_agoffset(method, l::Lookup, scale) = _agoffset(locus(l), scale)
_agoffset(locus::Start, scale) = 0
_agoffset(locus::End, scale) = scale - 1
_agoffset(locus::Center, scale) = scale ÷ 2

function _endpoints(method, l::Lookup, scale)
    start = firstindex(l) + _agoffset(method, l, scale)
    stop = (length(l) ÷ scale) * scale
    return start, stop
end
