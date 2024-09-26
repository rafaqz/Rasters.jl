
const DimOrDimTuple = Union{Dimension,Tuple{Vararg{<:Dimension}}}
const IntOrIntTuple = Union{Int,Tuple{Vararg{<:Int}}}

struct Ag end
struct DisAg end

const SKIPMISSING_KEYWORD = """
- `skipmissing`: if `true`, any `missingval` will be skipped during aggregation, so that
    only areas of all missing values will be aggregated to `missingval`. If `false`, any
    aggregated area containing a `missingval` will be assigned `missingval`.
"""
const METHOD_ARGUMENT = """
- `method`: a function such as `mean` or `sum` that can combine the
  value of multiple cells to generate the aggregated cell, or a [`Locus`]($DDlocusdocs)
  like `Start()` or `Center()` that species where to sample from in the interval.
"""
const SCALE_ARGUMENT = """
- `scale`: the aggregation factor, which can be an integer, a tuple of integers
  for each dimension, or any `Dimension`, `Selector` or `Int` combination you can
  usually use in `getindex`. Using a `Selector` will determine the scale by the
  distance from the start of the index in the `src` array.
"""

"""
    aggregate(method, object, scale; kw...)

Aggregate a `Raster`, or all arrays in a `RasterStack` or `RasterSeries`, by `scale` using
`method`.

# Arguments

$METHOD_ARGUMENT
- `object`: Object to aggregate, like `AbstractRasterSeries`, `AbstractStack`,
  `AbstractRaster` or `Dimension`.
$SCALE_ARGUMENT

When the aggregation `scale` of is larger than the array axis, the length of the axis is used.

# Keywords

$SKIPMISSING_KEYWORD
$FILENAME_KEYWORD
$SUFFIX_KEYWORD
$PROGRESS_KEYWORD
$THREADED_KEYWORD

# Example

```jldoctest
using Rasters, RasterDataSources, Statistics, Plots
using Rasters: Center
st = read(RasterStack(WorldClim{Climate}; month=1))
ag = aggregate(Center(), st, (Y(20), X(20)); skipmissing=true, progress=false)
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
    progress=true, threaded=false, kw...
)
    f(A) = aggregate(method, A, scale; progress=false, kw...)
    A1 = f(first(series))
    dest = similar(series, typeof(A1))
    dest[1] = A1
    _run(eachindex(series), threaded, progress, "Aggregating series...") do i
        dest[i] = f(series[i])
    end

    return dest
end
function aggregate(method, stack::AbstractRasterStack{K}, scale;
    keys=keys(stack), filename=nothing, suffix=keys, progress=true, threaded=false, kw...
) where K
    src = layers(stack)
    dst_vec = Vector{Raster}(undef, length(K)) # Intentionally type unstable
    _run(1:length(K), threaded, progress, "Aggregating stack...") do i
        dst_vec[i] = aggregate(method, src[i], scale; filename, suffix=suffix[i], kw...)
    end
    dst_tuple = ntuple(i -> dst_vec[i], Val{length(K)}())

    return DD.rebuild_from_arrays(stack, dst_tuple)
end
function aggregate(method, src::AbstractRaster, scale;
    suffix=nothing, filename=nothing, progress=true, kw...
)
    dst = alloc_ag(method, src, scale; filename, suffix, kw...)
    aggregate!(method, dst, src, scale; progress, kw...)
    return dst
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
    aggregate!(method, dst::AbstractRaster, src::AbstractRaster, scale; skipmissing=false)

Aggregate raster `src` to raster `dst` by `scale`, using `method`.

# Arguments

$METHOD_ARGUMENT
$SCALE_ARGUMENT

When the aggregation `scale` of is larger than the array axis, the length of the axis is used.

# Keywords

$SKIPMISSING_KEYWORD
$PROGRESS_KEYWORD

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
    return dst
end
# Function/functor methods
function aggregate!(f, dst::AbstractRaster, src, scale;
    skipmissingval=false, skipmissing=skipmissingval, progress=true
)
    intscale = _scale2int(Ag(), dims(src), scale)
    # len = prod(intscale)
    l = upsample.(map(firstindex, axes(dst)), intscale)
    u = upsample.(map(lastindex, axes(dst)), intscale)
    checkbounds(src, l...)
    checkbounds(src, u...)
    # If a disk array, cache the src so we don't read too many times
    src_parent = isdisk(src) ? DiskArrays.cache(parent(src)) : parent(src)
    @inbounds broadcast!(dst, CartesianIndices(dst)) do I
        upper = upsample.(Tuple(I), intscale)
        lower = upper .+ intscale .- 1
        I = map(:, upper, lower)
        block = isdisk(src) ? src_parent[I...] : Base.unsafe_view(src_parent, I...)
        if skipmissing
            _reduce_skip(f, block, missingval(src), dst)
        else
            _reduce_noskip(f, block, missingval(src), dst)
        end
    end
    return dst
end

"""
    disaggregate(object, scale; kw...)

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

$FILENAME_KEYWORD
$SUFFIX_KEYWORD
$PROGRESS_KEYWORD
$THREADED_KEYWORD

Note: currently it is faster to aggregate over memory-backed arrays.
Use [`read`](@ref) on `src` before use where required.

"""
function disaggregate end
disaggregate(_, x, scale) = disaggregate(x, scale) # legacy
function disaggregate(series::AbstractRasterSeries, scale;
    progress=true, threaded=false, kw...
)
    A1 = disaggregate(series[1], scale; progress=false, kw...)
    dst = similar(series, typeof(A1))
    dst[1] = A1
    _run(eachindex(series)[2:end], threaded, progress, "Disaggregating series...") do i
        dst[i] = disaggregate(series[i], scale; progress=false, kw...)
    end
    return dst
end
function disaggregate(stack::AbstractRasterStack{K}, scale;
    keys=keys(stack), suffix=keys, filename=nothing, progress=true, threaded=false
) where K
    dst_vec = Vector{Raster}(undef, length(K))
    ls = layers(stack)
    _run(1:length(K), threaded, progress, "Disaggregating stack...") do i
        dst_vec[i] = disaggregate(ls[i], scale; filename, suffix=suffix[i])
    end
    dst_tuple = ntuple(i -> dst_vec[i], Val{length(K)}())
    return DD.rebuild_from_arrays(stack, dst_tuple)
end
function disaggregate(src::AbstractRaster, scale;
    suffix=nothing, filename=nothing, kw...
)
    dst = alloc_disag(Center(), src, scale; filename, suffix, kw...)
    disaggregate!(dst, src, scale)
end
function disaggregate(dim::Dimension, scale)
    rebuild(dim, disaggregate(locus, lookup(dim), scale))
end
function disaggregate(lookup::Lookup, scale)
    intscale = _scale2int(DisAg(), lookup, scale)
    intscale == 1 && return lookup

    len = length(lookup) * intscale
    step_ = step(lookup) / intscale
    start = lookup[1] - _agoffset(Start(), intscale) * step_
    stop = start + (len - 1)  * step_
    index = LinRange(start, stop, len)
    if lookup isa AbstractSampled
        sp = disaggregate(locus, span(lookup), intscale)
        return rebuild(lookup; data=index, span=sp)
    else
        return rebuild(lookup; data=index)
    end
end

disaggregate(span::Span, scale) = span
disaggregate(span::Regular, scale) = Regular(val(span) / scale)

"""
    disaggregate!(dst::AbstractRaster, src::AbstractRaster, filename, scale)

Disaggregate array `src` to array `dst` by some scale.

- `scale`: the aggregation factor, which can be an integer, a tuple of integers
  for each dimension, or any `Dimension`, `Selector` or `Int` combination you can
  usually use in `getindex`. Using a `Selector` will determine the scale by the
  distance from the start of the index in the `src` array.

Note: currently it is faster to aggregate over memory-backed arrays.
Use [`read`](@ref) on `src` before use where required.
"""
function disaggregate!(dst::AbstractRaster, src, scale)
    disaggregate!((locus,), dst, src, scale)
end
function disaggregate!(dst::AbstractRaster, src, scale)
    intscale = _scale2int(DisAg(), dims(src), scale)
    broadcast!(dst, CartesianIndices(dst)) do I
        val = src[(downsample.(Tuple(I), intscale))...]
        val === missingval(src) ? missingval(dst) : val
    end
end

# Allocate an array of the correct size to aggregate `A` by `scale`
alloc_ag(method, A::AbstractRaster, scale; kw...) = alloc_ag((method,), A, scale; kw...)
function alloc_ag(method::Tuple, A::AbstractRaster, scale;
    filename=nokw, suffix=nokw, skipmissingval=false, skipmissing=false, progress=false
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
    filename=nokw, suffix=nokw
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
_scale2int(x, dims::DimTuple, scale::Tuple{<:Pair,Vararg{Pair}}) =
    _scale2int(x, dims, Dimensions.pairs2dims(scale...))
_scale2int(x, dims::DimTuple, scale::NamedTuple) = 
    _scale2int(x, dims, Dimensions.kw2dims(scale))
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

# Fallback iterator
@propagate_inbounds function _reduce_noskip(f, block, mv, dst)
    for x in block
        _ismissing(x, mv) && return _missingval_or_missing(dst)
    end
    return f(block)
end
# Specialised fast paths
@propagate_inbounds function _reduce_noskip(::typeof(count), block, mv, dst)
    Missings.nonmissingtype(eltype(block)) <: Bool || throw(ArgumentError("`count` can only reduce rasters of Bool"))
    return _reduce_noskip(sum, block, mv, dst)
end
@propagate_inbounds function _reduce_noskip(::typeof(sum), block, mv, dst)
    agg = zero(eltype(block))
    for x in block
        _ismissing(x, mv) && return _missingval_or_missing(dst)
        agg += x
    end
    return agg
end
@propagate_inbounds function _reduce_noskip(::typeof(DD.Statistics.mean), block, mv, dst)
    agg = zero(eltype(block))
    n = 0
    for x in block
        _ismissing(x, mv) && return _missingval_or_missing(dst)
        n += 1
        agg += x
    end
    return agg / n
end

# Fallback iterator
@propagate_inbounds function _reduce_skip(f, block, mv, dst)
    for x in block
        _ismissing(x, mv) || return f((x for x in block if x !== mv))
    end
    return _missingval_or_missing(dst)
end
# Specialised fast paths
@propagate_inbounds function _reduce_skip(::typeof(count), block, mv, dst)
    Missings.nonmissingtype(eltype(block)) <: Bool || throw(ArgumentError("`count` can only reduce rasters of Bool"))
    return _reduce_skip(sum, block, mv, dst)
end
@propagate_inbounds function _reduce_skip(::typeof(sum), block, mv, dst)
    agg = zero(eltype(block))
    found = false
    for x in block
        _ismissing(x, mv) && continue
        found = true
        agg += x
    end
    return found ? agg : _missingval_or_missing(dst)
end
@propagate_inbounds function _reduce_skip(::typeof(DD.Statistics.mean), block, mv, dst)
    agg = zero(eltype(block))
    found = false
    n = 0
    for x in block
        _ismissing(x, mv) && continue
        found = true
        n += 1
        agg += x
    end
    return found ? agg / n : _missingval_or_missing(dst)
end

_ismissing(x, mv) = ismissing(x) || x === mv 
