
const DimOrDimTuple = Union{Dimension,Tuple{Vararg{Dimension}}}
const IntOrIntTuple = Union{Int,Tuple{Vararg{Int}}}

struct Ag end
struct DisAg end

const SKIPMISSING_KEYWORD = """
- `skipmissing`: if `true`, any `missingval` will be skipped during aggregation, so that
    only areas of all missing values will be aggregated to `missingval(dst)`. If `false`,
    aggregated areas containing one or more `missingval` will be assigned `missingval`. 
    `false` by default. `skipmissing` behaviour is independent of function `f`, which is
    only applied to completely non-missing values.
"""
const METHOD_ARGUMENT = """
- `method`: a function such as `mean` or `sum` that can combine the
  value of multiple cells to generate the aggregated cell, or a [`Locus`]($DDlocusdocs)
  like `Start()` or `Center()` that species where to sample from in the interval.
"""
const SCALE_ARGUMENT = """
- `scale`: the aggregation factor, which can be an `Int`, a `Tuple` of `Int`
  for each dimension, or a `:` colon to mean the whole dimension. 
  You can also use any `Dimension`, `Selector` or `Int` combination you can
  usually use in `getindex`. `Tuple` of `Pair` or `NamedTuple` where keys are dimension names
  will also work. Using a `Selector` will determine the scale by the distance from the start 
  of the index. Selectors will find the first offset and repeat the same aggregation size for the rest.
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
$VERBOSE_KEYWORD

# Example

```jldoctest
using Rasters, RasterDataSources, Statistics, Plots
import ArchGDAL
using Rasters: Center
st = RasterStack(WorldClim{Climate}; month=1)
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
function aggregate(method, series::AbstractRasterSeries, scale;
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
    return alloc_ag(method, src, scale; filename, suffix, kw...) do dst
        aggregate!(method, dst, src, scale; progress, kw...)
    end
end
aggregate(method, d::Dimension, scale) = rebuild(d, aggregate(method, lookup(d), scale))
aggregate(method, l::Lookup, scale::Colon) = aggregate(method, l, length(l)) 
aggregate(method, l::Lookup, scale::Nothing) = aggregate(method, l, 1) 
function aggregate(method, l::Lookup, scale::Int)
    intscale = _scale2int(Ag(), l, scale)
    if issampled(l) && isordered(l) && isregular(l)
        start, stop = _endpoints(method, l, intscale)
        sp = aggregate(span(l), scale)
        return rebuild(l; data = start:val(sp):stop, span=sp)
    else
        # Categorical and Unordered lookups are just broken 
        # by aggregate, so use NoLookup
        return NoLookup(Base.OneTo(length(l) ÷ intscale))
    end
end
aggregate(span::Span, scale) = span
aggregate(span::Regular, scale) = Regular(val(span) * scale)

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
$VERBOSE_KEYWORD

Note: currently it is _much_ faster to aggregate over memory-backed source
arrays. Use [`read`](@ref) on `src` before use where required.
"""
function aggregate!(locus::Locus, dst::AbstractRaster, src, scale; kw...)
    aggregate!((locus,), dst, src, scale)
end
function aggregate!(loci::Tuple{Locus,Vararg}, dst::AbstractRaster, src, scale; 
    verbose=true, kw...
)
    comparedims(dst, src; length=false)
    intscale = _scale2int(Ag(), dims(src), scale; verbose)
    offsets = ceil.(Int, _agoffset.(loci, (ForwardOrdered(),), intscale))
    # Cache the source if its a disk array
    src1 = isdisk(src) ? DiskArrays.cache(src) : src
    # Broadcast will make the dest arrays chunks when needed
    broadcast!(dst, CartesianIndices(dst)) do I
        # Upsample for each pixel. Possibly slighly inneficient
        # for large aggregations but its simple
        J = upsample.(Tuple(I), intscale) .+ offsets
        val = src1[J...]
        _ismissing(val, missingval(src)) ? missingval(dst) : val
    end
    return dst
end
# Function/functor methods
function aggregate!(f, dst::AbstractRaster, src, scale;
    skipmissingval=false, skipmissing=skipmissingval, progress=true, verbose=true
)
    @show dims(dst)
    comparedims(dst, src; length=false)
    all(Lookups.isaligned, lookup(src)) || 
        throw(ArgumentError("Currently only grid-alligned dimensions can be aggregated. Make a Rasters.jl Github issue if you need to aggregate with transformed dims"))

    intscale = _scale2int(Ag(), dims(src), scale; verbose)
    l = upsample.(map(firstindex, axes(dst)), intscale)
    u = upsample.(map(lastindex, axes(dst)), intscale)
    checkbounds(src, l...)
    checkbounds(src, u...)
    # If a disk array, cache the src so we don't read too many times
    # where src and dest chunks really don't align this may not be efficient
    open(src) do src1
        src2 = isdisk(src1) ? DiskArrays.cache(src1) : src1
        # Broadcast will make the dest arrays chunks when needed
        @inbounds broadcast!(dst, CartesianIndices(dst)) do I
            # Upsample the lower bounds of the source range
            lower = upsample.(Tuple(I), intscale)
            upper = map(lower, intscale) do u, i
                isnothing(i) ? u : u + i - 1
            end
            I = map(:, lower, upper)
            # Read a block from the source array to reduce
            block = isdisk(src2) ? parent(src2)[I...] : Base.view(parent(src2), I...)
            # Reduce with or without skipmissing
            if skipmissing
                # With skipmissing `f` will recieve an iterator
                _reduce_skip(f, block, missingval(src), dst)
            else
                # Without skipmissing `f` will recieve a Raster
                _reduce_noskip(f, block, missingval(src), dst)
            end
        end
    end
    return dst
end

"""
    disaggregate(object, scale; kw...)

Disaggregate array, or all arrays in a stack or series, by some scale.

# Arguments

- `object`: Object to aggregate, like `AbstractRasterSeries`, `AbstractStack`,
  `AbstractRaster`, `Dimension` or `Lookup`.
$SCALE_ARGUMENT

# Keywords

$FILENAME_KEYWORD
$SUFFIX_KEYWORD
$PROGRESS_KEYWORD
$THREADED_KEYWORD
- `lazy`: A `Bool` specifying if to disaggregate lazily. Defaults to `false`

Note: currently it is _much_ faster to disaggregate over a memory-backed 
source array. Use [`read`](@ref) on `src` before use where required.
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
    keys=keys(stack), suffix=keys, progress=true, threaded=false, kw...
) where K
    dst_vec = Vector{Raster}(undef, length(K))
    ls = layers(stack)
    _run(1:length(K), threaded, progress, "Disaggregating stack...") do i
        dst_vec[i] = disaggregate(ls[i], scale; suffix=suffix[i], kw...)
    end
    dst_tuple = ntuple(i -> dst_vec[i], Val{length(K)}())
    return DD.rebuild_from_arrays(stack, dst_tuple)
end
function disaggregate(src::AbstractRaster, scale;
    suffix=nothing, filename=nothing, lazy = false, kw...
)
    if lazy
        return view_disaggregate(src, scale)
    else
        return alloc_disag(Center(), src, scale; filename, suffix, kw...) do dst
            disaggregate!(dst, src, scale)
        end
    end
end
function disaggregate(dim::Dimension, scale)
    rebuild(dim, disaggregate(locus, lookup(dim), scale))
end
function disaggregate(l::Lookup, scale)
    intscale = _scale2int(DisAg(), l, scale)
    intscale == 1 && return l

    len = length(l) * intscale
    step_ = step(l) / intscale
    start = first(l) - _agoffset(l, intscale) * step_
    stop = start + (len - 1)  * step_
    index = LinRange(start, stop, len)
    if l isa AbstractSampled
        sp = disaggregate(locus, span(l), intscale)
        rebuild(l; data=index, span=sp)
    else
        rebuild(l; data=index)
    end
end

disaggregate(span::Span, scale) = span
disaggregate(span::Regular, scale) = Regular(val(span) / scale)

function view_disaggregate(A, scale)
    intscale = _scale2int(DisAg(), dims(A), scale)
    dims_ = disaggregate.((Center(),), dims(A), intscale)
    indices = map((a, i) -> repeat(a; inner =i), axes(A), intscale)
    rebuild(A; data = view(parent(A), indices...), dims = dims_)
end
"""
    disaggregate!(dst::AbstractRaster, src::AbstractRaster, scale)

Disaggregate array `src` to array `dst` by some scale.

$SCALE_ARGUMENT
"""
function disaggregate!(dst::AbstractRaster, src, scale)
    intscale = _scale2int(DisAg(), dims(src), scale)
    # For now we just read a DiskArray
    open(src) do src1
        src2 = isdisk(src) ? DiskArrays.cache(src1) : src1
        # Fast path for Array backed data
        for I in CartesianIndices(src2)
            lower = upsample.(Tuple(I), intscale)
            upper = map(lower, intscale, size(dst)) do l, is, s
                min(l + is - 1, s) 
            end
            ranges = map(:, lower, upper)
            val = src2[I]
            val1 = _ismissing(src, val) ? missingval(dst) : val
            dst[ranges...] .= (val1,)
        end
    end
    return dst
end

# Allocate an array of the correct size to aggregate `A` by `scale`
alloc_ag(f, method, A::AbstractRaster, scale; kw...) = alloc_ag(f, (method,), A, scale; kw...)
function alloc_ag(f, method::Tuple, A::AbstractRaster, scale;
    filename=nokw, suffix=nokw, skipmissingval=false, skipmissing=false, progress=false, verbose=false
)
    intscale = _scale2int(Ag(), dims(A), scale; verbose=false)
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
    return create(f, filename, T, dims_; name=name(A), suffix, missingval=mv)
end

# Allocate an array of the correct size to disaggregate `A` by `scale`
function alloc_disag(f, method, A::AbstractRaster, scale; kw...)
    alloc_disag(f, (method,), A, scale; kw...)
end
function alloc_disag(f, method::Tuple, A::AbstractRaster, scale;
    filename=nokw, suffix=nokw
)
    intscale = _scale2int(DisAg(), dims(A), scale; verbose=false)
    dims_ = map(dims(A), intscale) do d, i
        disaggregate(method, d, i)
    end
    # Dim aggregation determines the array size
    sze = map(length, dims_)
    T = ag_eltype(method, A)
    mv = missingval(A) isa Nothing ? nothing : convert(T, missingval(A))
    return create(f, filename, T, dims_; name=name(A), suffix, missingval=mv)
end

# Handle how methods like `mean` can change the type
ag_eltype(method::Tuple{<:Locus,Vararg}, A) = eltype(A)
function ag_eltype(method::Tuple{<:Any}, A)
    method_returntype = typeof(method[1](zero(eltype(A))))
    promote_type(eltype(A), method_returntype)
end

# Convert indices from the aggregated array to the larger original array.
upsample(index::Int, scale::Int) = (index - 1) * scale + 1
upsample(index::Int, scale::Nothing) = index

# Convert indices from the original array to the aggregated array.
downsample(index::Int, scale::Int) = (index - 1) ÷ scale + 1
downsample(index::Int, scale::Nothing) = index

# Convert scale or tuple of scale to integer using dims2indices
@inline function _scale2int(x, dims::DimTuple, scale::DimTuple; verbose=true)
    map(dims, DD.sortdims(scale, dims)) do d, s
        if isnothing(s) 
            1
        else
            i = dims2indices(d, s)
            # Swap Colon as all to Colon as 1 (1 is the no-change option here)
            s = i isa Colon ? length(d) : i
            _scale2int(x, d, s)
        end
    end
end
@inline function _scale2int(x, dims::DimTuple, scale::Tuple; verbose=true)
    map(dims, DD.dims2indices(dims, scale)) do d, s
        _scale2int(x, d, s)
    end
end
@inline _scale2int(x, dims::DimTuple, scale::Tuple{<:Pair,Vararg{Pair}}; verbose=true) =
    _scale2int(x, dims, Dimensions.pairs2dims(scale...); verbose)
@inline _scale2int(x, dims::DimTuple, scale::NamedTuple; verbose=true) = 
    _scale2int(x, dims, Dimensions.kw2dims(scale); verbose)
@inline _scale2int(x, dims::DimTuple, scale::Dimension; verbose=true) = 
    _scale2int(x, dims, (scale,); verbose)
@inline function _scale2int(x, dims::DimTuple, scale::Int; verbose=true) 
    # If there are other dimensions, we skip categorical dims
    vals = map(dims) do d
        if iscategorical(d) || !isordered(d) 
            name(d), nothing
        else
            name(d), _scale2int(x, d, scale)
        end
    end
    nskipped = count(isnothing ∘ last, vals)
    if nskipped == length(dims)
        example = join(map(((n, d),) -> "$(name(d))=$(n + 1),", enumerate(dims)), ' ')
        # If all dims are categorical we error
        throw(ArgumentError("All dimensions are Categorical. To aggregate anyway, list scale explicity for each dimension, e.g. ($example)"))
    end
    if verbose && nskipped > 0
        scaleddims = join((v[1] for v in vals if v[2] isa Int), ", ", " and ")
        skippeddims = join((v[1] for v in vals if isnothing(v[2])), ", ", " and ")
        @info """
            Aggregating $scaleddims by $scale. $(skippeddims == "" ? "" : skippeddims) 
            skipped due to being `Categorical` or `Unordered`. 
            Specify all scales explicitly in a Tuple or NamedTuple to aggregate these anyway.  
            """
    end
    return map(last, vals)
end
@inline _scale2int(x, dims::DimTuple, scale::Colon; verbose=true) = 
    _scale2int(x, dims, map(_ -> Colon(), dims)) 
@inline _scale2int(x, d, scale::Colon) = length(d)
@inline _scale2int(x, dim::Dimension, scale::Int) = _scale2int(x, lookup(dim), scale)
@inline _scale2int(::Ag, l::Lookup, scale::Int) = scale > length(l) ? length(l) : scale
@inline _scale2int(::DisAg, l::Lookup, scale::Int) = scale

_agoffset(locus::Locus, l::Lookup, scale::Int) = _agoffset(locus, scale)
_agoffset(method, l::Lookup, scale::Int) = _agoffset(l, scale)
_agoffset(l::Lookup, scale::Int) = _agoffset(locus(l), order(l), scale)
_agoffset(x, scale::Colon) = 0
_agoffset(locus::Start, ::ForwardOrdered, scale::Int) = 0
_agoffset(locus::End, ::ForwardOrdered, scale::Int) = scale - 1
_agoffset(locus::Start, ::ReverseOrdered, scale::Int) = scale - 1
_agoffset(locus::End, ::ReverseOrdered, scale::Int) = 0
_agoffset(locus::Center, ::Ordered, scale::Int) = (scale-1)/2

_endpoints(_, l::Lookup, scale::Int) = _endpoints(locus(l), l, scale)
function _endpoints(locus::Locus, l::Lookup, scale::Int)
    offset = step(l)*_agoffset(locus, order(l), scale)
    start = first(l) + offset
    stop = l[(length(l) ÷ scale)*scale] + offset
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
