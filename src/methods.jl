const GeoStackOrArray = Union{AbstractGeoStack,AbstractGeoArray}
const GeoSeriesOrStack = Union{AbstractGeoSeries,AbstractGeoStack}

"""
    replace_missing(a::AbstractGeoArray, newmissingval)
    replace_missing(a::AbstractGeoStack, newmissingval)

Replace missing values in the array or stack with a new missing value,
also updating the `missingval` field/s.

A `GeoArray` containing a newly allocated `Array` is always returned,
even when the missing value matches the current value.
"""
replace_missing(x; missingval=missing) = replace_missing(x, missingval)
function replace_missing(A::AbstractGeoArray, missingval=missing)
    newdata = if ismissing(GeoData.missingval(A))
        if ismissing(missingval)
            copy(parent(read(A)))
        else
            collect(Missings.replace(parent(A), missingval))
        end
    else
        replace(parent(A), GeoData.missingval(A) => missingval)
    end
    rebuild(A; data=newdata, missingval=missingval)
end
replace_missing(x::GeoSeriesOrStack, args...) = map(A -> replace_missing(A, args...), x)

"""
    boolmask(A::AbstractArray, [missingval])

Create a mask array of `Bool` values, from any AbstractArray. For `AbstractGeoArray`
the default `missingval` is `missingval(A)`, for all other `AbstractArray`s
it is `missing`.

The array returned from calling `boolmask` on a `AbstractGeoArray` is a
[`GeoArray`](@ref) with the same size and fields as the oridingl array
"""
function boolmask end
boolmask(A::AbstractGeoArray) =
    rebuild(A; data=boolmask(A, missingval(A)), missingval=false, name=:Bool_mask)
boolmask(A::AbstractArray, missingval::Missing=missing) = (a -> !ismissing(a)).(parent(A))
boolmask(A::AbstractArray, missingval) =
    if isnan(missingval)
        (a -> !isnan(a)).(parent(A))
    else
        (a -> !isapprox(a, missingval)).(parent(A))
    end

"""
    missingmask(A::AbstractArray, [missingval])

Create a mask array of `missing` or `true` values, from any AbstractArray.
For `AbstractGeoArray` the default `missingval` is `missingval(A)`,
for all other `AbstractArray`s it is `missing`.

The array returned from calling `boolmask` on a `AbstractGeoArray` is a
[`GeoArray`](@ref) with the same size and fields as the oridingl array
"""
function missingmask end
missingmask(A::AbstractGeoArray) =
    rebuild(A; data=missingmask(A, missingval(A)), missingval=missing, name=:Missing_mask)
missingmask(A::AbstractArray, missingval::Missing=missing) =
    (a -> ismissing(a) ? missing : true).(parent(A))
missingmask(A::AbstractArray, missingval) =
    if isnan(missingval)
        (a -> isnan(a) ? missing : true).(parent(A))
    else
        (a -> isapprox(a, missingval) ? missing : true).(parent(A))
    end

"""
    mask(A::AbstractGeoArray; to; missingval=missingval(A))

Return a new array with values of `A` masked by the missing values of `to`,
optionally setting a new `misssingval`.

$EXPERIMENTAL
"""
function mask end
function mask(A::AbstractGeoArray; to, missingval=missingval(A))
    missingval = missingval isa Nothing ? missing : missingval
    return mask!(replace_missing(A, missingval); to, missingval)
end
mask(xs::GeoSeriesOrStack; kw...) = map(x -> mask(x; kw...),  xs)
"""
    mask(A::AbstractGeoArray, polygon; order, [, filename])`

Masks the array `A` by a polygon. 

Creates a new array, where points falling outside the polygon have 
been replaced by `missingval(A)`.

$EXPERIMENTAL
"""
function mask(A::AbstractGeoArray, poly::AbstractVector; order=(XDim, YDim), filename=nothing)
    broadcast(DimKeys(A)) do ds
        _polymask(src, dst, ds)
    end
end

"""
    mask!(A; to::AbstractArray, missingval=missing)

Mask `A` by the missing values of `to`, optionally setting a new `misssingval`.

$EXPERIMENTAL
"""
function mask!(A::AbstractGeoArray; to, missingval=missingval(A))
    missingval isa Nothing && throw(ArgumentError("Array has no `missingval`. Pass a `missingval` keyword compatible with the type, or use `rebuild(A; missingval=somemissingval)` to set it."))
    dimwise!(A, A, to) do a, t
        t === GeoData.missingval(to) ? missingval : a
    end
    return A
end
function mask!(xs::GeoSeriesOrStack; kw...)
    map(x -> mask!(x; kw...),  xs)
    return xs
end
function mask!(dst, src::AbstractGeoArray, poly::AbstractVector; order=(XDim, YDim), filename=nothing)
    broadcast!(dst, DimKeys(src)) do ds
        _polymask(dst, src, ds)
    end
end

function _polymask(dst, src, dims)
    x = src[dims...]
    if x === missingval(src)
        missingval(dst)
    elseif isinside(map(val, dims), poly)
        x
    else
        missingval(dst)
    end
end

const Pt{T<:Real} = NTuple{N,T} where N

isinside(sel::At, poly::AbstractVector) = isinside(val(sel), poly)
function isinside(r::Pt, poly::AbstractVector)
    # An implementation of Hormann-Agathos (2001) Point in Polygon algorithm
    # See: http://www.sciencedirect.com/science/article/pii/S0925772101000128
    # Code segment adapted from PolygonClipping.jl
    c = false
    detq(q1,q2,r) = (q1[1] - r[1]) * (q2[2] - r[2]) - (q2[1] - r[1]) * (q1[2] - r[2])

    for i in eachindex(poly)[2:end]
        q2 = poly[i]
        q1 = poly[i - 1]
        if q1 == r
            @warn("point on polygon vertex - returning false")
            return false
        end
        if q2[2] == r[2]
            if q2[1] == r[1]
                @warn("point on polygon vertex - returning false")
                return false
            elseif (q1[2] == r[2]) && ((q2[1] > x) == (q1[1] < r[1]))
                @warn("point on edge - returning false")
                return false
            end
        end
        if (q1[2] < r[2]) != (q2[2] < r[2]) # crossing
            if q1[1] >= r[1]
                if q2[1] > r[1]
                    c = !c
                elseif ((detq(q1,q2,r) > 0) == (q2[2] > q1[2])) # right crossing
                    c = !c
                end
            elseif q2[1] > r[1]
                if ((detq(q1,q2,r) > 0) == (q2[2] > q1[2])) # right crossing
                    c = !c
                end
            end
        end
    end
    return c
end

"""
    classify(x, pairs; lower, upper, others)
    classify(x, pairs...; lower, upper, others)

Create a new array with values in `x` classified by the values in `pairs`, where each
pair contains a value and a replacement, a tuple of lower and upper range and a replacement,
or a Tuple of `Fix2` like `>(x)`. If `Fix2` is not used, the `lower` and `upper` keywords
determined which comparison operator is used for the lower and upper tuple values
- e.g. `<` or `<=` and `>` or `>=`.

If `others` is set other values not covered in `pairs` will be set to that values.

$EXPERIMENTAL
"""
function classify end
classify(A::AbstractGeoArray, pairs::Pair...; kw...) = classify(A, pairs; kw...)
function classify(A::AbstractGeoArray, pairs; lower=(>=), upper=(<), others=nothing)
    broadcast(A) do x
        _classify(x, pairs, lower, upper, others, missingval(A))
    end
end
classify(xs::GeoSeriesOrStack, values; kw...) = map(x -> classify(x, values; kw...),  xs)

"""
    classify!(x, pairs...; lower, upper, others)
    classify!(x, pairs; lower, upper, others)
    classify!(x, pairs; lower, upper, others)

Classify values in `x` by the values in `pairs`, where each pair contains a value
and a replacement, a tuple of lower and upper range and a replacement, or a Tuple
of `Fix2` like `>(x)`. If `Fix2` is not used, the `lower` and `upper` keywords
determined which comparison operator is used for the lower and upper tuple values
- e.g. `<` or `<=` and `>` or `>=`.

If `others` is set other values not covered in `pairs` will be set to that values.

$EXPERIMENTAL
"""
classify!(A::AbstractGeoArray, pairs::Pair...; kw...) = classify!(A, pairs; kw...)
function classify!(A::AbstractGeoArray, pairs; lower=(>=), upper=(<), others=nothing)
    broadcast!(A, A) do x
        _classify(x, pairs, lower, upper, others, missingval(A))
    end
end
function classify!(xs::GeoSeriesOrStack; kw...)
    map(x -> classify!(x; kw...),  xs)
    return xs
end

# _classify
# Classify single values
function _classify(x, pairs, lower, upper, others, missingval)
    x === missingval && return x
    # Use a fold instead of a loop, for type stability
    found = foldl(pairs; init=nothing) do found, (find, replace)
        if found isa Nothing && _compare(find, x, lower, upper)
            replace
        else
            found
        end
    end
    if found isa Nothing
        if others isa Nothing
            return x
        else
            return others
        end
    else
        return found
    end
end
function _classify(x, pairs::AbstractMatrix, lower, upper, others, missingval)
    x === missingval && return x
    found = false
    if size(pairs, 2) == 2
        for i in 1:size(pairs, 1)
            find = pairs[i, 1]
            if _compare(find, x, lower, upper)
                x = pairs[i, 2]
                found = true
                break
            end
        end
    elseif size(pairs, 2) == 3
        for i in 1:size(pairs, 1)
            find = pairs[i, 1], pairs[i, 2]
            if _compare(find, x, lower, upper)
                x = pairs[i, 3]
                found = true
                break
            end
        end
    else
        throw(ArgumentError("pairs Array must be a N*2 or N*3 matrix"))
    end
    if !found && !(others isa Nothing)
        x = others
    end
    return x
end

_compare(find, x, lower, upper) = find === x
_compare(find::Base.Fix2, x, lower, upper) = find(x)
_compare((l, u)::Tuple, x, lower, upper) = lower(x, l) && upper(x, u)
_compare((l, u)::Tuple{<:Base.Fix2,<:Base.Fix2}, x, lower, upper) = l(x) && u(x)


"""
    crop(layers::AbstractGeoArray...)
    crop(layers::Union{NamedTuple,Tuple})
    crop(A::Union{AbstractGeoArray,AbstractGeoStack}; to::Tuple)

Crop multiple [`AbstractGeoArray`](@ref) to match the size
of the smallest one for any dimensions that are shared.

Otherwise crop to the size of the keyword argument `to`. This can be a
`Tuple` of `Dimension` or any object that will return one from `dims(to)`.

$EXPERIMENTAL
"""
function crop end
crop(l1, l2, ls::AbstractGeoArray...; kw...) = crop((l1, l2, ls); kw...)
function crop(layers::Union{Tuple,NamedTuple}; to=_smallestdims(layers))
    map(l -> crop(l; to), layers)
end
crop(x::GeoStackOrArray; to) = _crop_to(x, to)

# crop `A` to values of dims of `to`
_crop_to(A::GeoStackOrArray, to) = _crop_to(A, dims(to))
function _crop_to(x::GeoStackOrArray, to::Tuple)
    # Create selectors for each dimension
    # `Between` the bounds of the dimension
    selectors = map(to) do d
        DD.basetypeof(d)(Between(DD.bounds(d)))
    end
    # Take a view of the `Between` selectors
    return _without_mapped_crs(x) do a
        view(a, selectors...)
    end
end

# Get the smallest dimensions in a tuple of AbstractGeoArray
function _smallestdims(layers)
    # Combine the dimensions of all layers
    dims = DD.combinedims(layers...; check=false)
    # Search through all the dimensions choosing the shortest
    alldims = map(DD.dims, layers)
    return map(dims) do d
        matchingdims = map(ds -> DD.dims(ds, (d,)), alldims)
        reduce(matchingdims) do a, b
            _choose(_shortest, a, b)
        end |> first
    end
end

"""
    extend(layers::AbstractGeoArray...)
    extend(layers::Union{NamedTuple,Tuple})
    extend(A::Union{AbstractGeoArray,AbstractGeoStack}; to)

Extend multiple [`AbstractGeoArray`](@ref) to match the area covered by all.
A single `AbstractGeoArray` can be extended by passing the new `dims` tuple
as the second argument.

$EXPERIMENTAL
"""
function extend end
function extend(l1::AbstractGeoArray, l2::AbstractGeoArray, ls::AbstractGeoArray...; kw...)
    extend((l1, l2, ls...); kw...)
end
function extend(layers::Union{NamedTuple,Tuple}; to=_largestdims(layers))
    # Extend all layers to `to`, by default the _largestdims
    map(l -> extend(l; to), layers)
end
extend(x::GeoStackOrArray; to=dims(x)) = _extend_to(x, to)

_extend_to(x::GeoStackOrArray, to) = _extend_to(x, dims(to))
function _extend_to(A::AbstractGeoArray, to::Tuple)
    sze = map(length, to)
    T = eltype(A)
    # Create a new extended array
    newdata = similar(parent(A), T, sze)
    # Fill it with missing/nodata values
    newdata .= missingval(A)
    # Rebuild the original object with larger data and dims.
    newA = rebuild(A; data=newdata, dims=to)
    # Calculate the range of the old array in the extended array
    ranges = map(dims(A), to) do d, nd
        start = DD.sel2indices(nd, Near(first(d)))
        stop = DD.sel2indices(nd, Near(last(d)))
        start <= stop ? (start:stop) : (stop:start)
    end
    # Copy the original data to the new array
    copyto!(
        parent(newA), CartesianIndices((ranges...,)),
        parent(read(A)), CartesianIndices(A)
    )
    return newA
end
_extend_to(st::AbstractGeoStack, to::Tuple) = map(A -> _extend_to(A, to), st)

# Get the largest dimensions in a tuple of AbstractGeoArray
function _largestdims(layers)
    dims = DD.combinedims(layers...; check=false)
    alldims = map(DD.dims, layers)
    return map(dims) do d
        matchingdims = map(ds -> DD.dims(ds, (d,)), alldims)
        reduce(matchingdims) do a, b
            _choose(_longest, a, b)
        end |> first
    end
end

# Choose a dimension from either missing dimension
# (empty Tuple) or a comparison between two 1-Tuples
_choose(f, ::Tuple{}, ::Tuple{}) = ()
_choose(f, ::Tuple{}, (b,)::Tuple) = (b,)
_choose(f, (a,)::Tuple, ::Tuple{}) = (a,)
_choose(f, (a,)::Tuple, (b,)::Tuple) = (f(a, b) ? a : b,)

# Choose the shortest or longest dimension
_shortest(a, b) = length(a) <= length(b)
_longest(a, b) = length(a) >= length(b)

"""
    trim(A::AbstractGeoArray; dims::Tuple, pad::Int)

Trim `missingval` from `A` for axes in dims, returning a view of `A`.

By default `dims=(X, Y)`, so trimming keeps the area of `X` and `Y`
that contains non-missing values along all other dimensions.

The trimmed size will be padded by `pad` on all sides, although
padding will not be added beyond the original extent of the array.

$EXPERIMENTAL
"""
function trim(A::GeoStackOrArray; dims::Tuple=(X(), Y()), pad::Int=0)
    # Get the actual dimensions in their order in the array
    dims = commondims(A, dims)
    # Get the range of non-missing values for each dimension
    ranges = _trimranges(A, dims)
    # Add paddding
    padded = map(ranges, map(d -> size(A, d), dims)) do r, l
        max(first(r)-pad, 1):min(last(r)+pad, l)
    end
    dims = map(rebuild, dims, padded)
    return view(A, dims...)
end

# Tracks the status of an index for some subset of dimensions of an Array
# This lets us track e.g. the X/Y indices that have only missing values
# accross all other dimensions.
# This is a hack to work with DiskArrays broadcast chunking without allocations.
struct AxisTrackers{N,Tr,D,TD} <: AbstractArray{Bool,N}
    tracking::Tr
    dims::D
    trackeddims::TD
end
function AxisTrackers(tracking::T, dims::D, trackeddims::TD) where {T,D,TD}
    AxisTrackers{length(dims),T,D,TD}(tracking, dims, trackeddims)
end
function AxisTrackers(dims::Tuple, trackeddims::Tuple)
    tracking = map(trackeddims) do td
        (_ -> false).(td)
    end
    return AxisTrackers(tracking, dims, trackeddims)
end

Base.axes(A::AxisTrackers) = map(d -> axes(d, 1), A.dims)
Base.size(A::AxisTrackers) = map(length, A.dims)
Base.getindex(A::AxisTrackers, I...) = map(getindex, A.tracking, _trackedinds(I)) |> any
function Base.setindex!(A::AxisTrackers, x, I::Int...)
    map(A.tracking, _trackedinds(A, I)) do axis, i
        axis[i] |= x
    end
end

function _trackedinds(A, I)
    # Wrap indices in dimensions so we can sort and filter them
    Id = map((d, i) -> DD.basetypeof(d)(i), A.dims, I)
    # Get just the tracked dimensions
    Itracked = dims(Id, A.trackeddims)
    # Get the indices for the tracked dimensions
    return map(val, Itracked)
end

# Get the ranges to trim to for dimensions in `dims`
function _trimranges(A, targetdims)
    # Broadcast over the array and tracker to mark axis indices
    # as being missing or not
    trackers = AxisTrackers(dims(A), targetdims)
    _update!(trackers, A)
    # Get the ranges that contain all non-missing values
    cropranges = map(trackers.tracking) do a
        f = findfirst(a)
        l = findlast(a)
        f = f === nothing ? firstindex(a) : f
        l = l === nothing ? lastindex(a) : l
        f:l
    end
    return cropranges
end

_update!(tr::AxisTrackers, A::AbstractGeoArray) = tr .= A .!== missingval(A)
_update!(tr::AxisTrackers, st::AbstractGeoStack) = map(A -> tr .= A .!== missingval(A), st)

"""
    slice(A::Union{AbstractGeoArray,AbstractGeoStack,AbstracGeoSeries}, dims)

Slice an object along some dimension/s, lazily using `view`.

For a single `GeoArray` or `GeoStack` this will return a `GeoSeries` of
`GeoArray` or `GeoStack` that are slices along the specified dimensions.

For a `GeoSeries`, the output is another series where the child objects are sliced and the
series dimensions index is now of the child dimensions combined. `slice` on a `GeoSeries`
with no dimensions will slice along the dimensions shared by both the series and child object.

$EXPERIMENTAL
"""
slice(x::GeoStackOrArray, dims) = slice(x, (dims,))
# Slice an array or stack into a series
function slice(x::GeoStackOrArray, dims::Tuple)
    # Make sure all dimensions in `dims` are in `x`
    all(hasdim(x, dims)) || _errordimsnotfound(dims, DD.dims(x))
    # Define dimensions and data for the sliced GeoSeries
    seriesdims = DD.dims(x, dims)
    # series data is a generator of view slices
    seriesdata = map(DimIndices(seriesdims)) do ds
        view(x, ds...)
    end
    return GeoSeries(seriesdata, seriesdims)
end
# Slice an existing series into smaller slices
slice(ser::AbstractGeoSeries, dims) = cat(map(x -> slice(x, dims), ser)...; dims=dims)

@noinline _errordimsnotfound(targets, dims) =
    throw(ArgumentError("Dimensions $(map(DD.dim2key, targets)) were not found in $(map(DD.dim2key, dims))"))

# By default, combine all the GeoSeries dimensions and return a GeoArray or GeoStack
combine(ser::AbstractGeoSeries) = combine(ser, dims(ser))
# Fold over all the dimensions, combining the series one dimension at a time
combine(ser::AbstractGeoSeries, dims::Tuple) = foldl(combine, dims; init=ser)
# Slice the N-dimensional series into an array of 1-dimensional
# series, and combine them, returning a new series with 1 less dimension.
function combine(ser::AbstractGeoSeries{<:Any,M}, dim::Union{Dimension,DD.DimType,Val,Symbol}) where M
    od = otherdims(ser, dim)
    slices = map(d -> view(ser, d...), DimIndices(od))
    newchilren = map(s -> combine(s, dim), slices)
    return rebuild(ser; data=newchilren, dims=od)
end
# Actually combine a 1-dimensional series with `cat`
function combine(ser::AbstractGeoSeries{<:Any,1}, dim::Union{Dimension,DD.DimType,Val,Symbol})
    dim = DD.dims(ser, dim)
    D = DD.basetypeof(dim)
    x = foldl(ser) do acc, x
        # May need to reshape to match acc
        cat(acc, _maybereshape(x, acc, dim); dims=D)
    end
    return set(x, D => dims(ser, dim))
end

function _maybereshape(A::AbstractGeoArray{<:Any,N}, acc, dim) where N
    if ndims(acc) != ndims(A)
        newdata = reshape(parent(A), Val{N+1}())
        d = if hasdim(refdims(A), dim)
            dims(refdims(A), dim)
        else
            DD.basetypeof(dim)(1:1; mode=NoIndex())
        end
        newdims = (DD.dims(A)..., d)
        return rebuild(A; data=newdata, dims=newdims)
    else
        return A
    end
end
function _maybereshape(st::AbstractGeoStack, acc, dim)
    map((s, a) -> _maybereshape(s, a, dim), st, acc)
end

"""
    chunk(A::AbstractGeoArray)

Creat a GeoSeries of arrays matching the chunks of a chunked array.

This may be useful for parallel or larger than memory applications.

$EXPERIMENTAL
"""
function chunk(A::AbstractGeoArray)
    # Get the index of each chunk of A
    gc = DiskArrays.eachchunk(A)
    ci = CartesianIndices(gc.chunkgridsize)
    # Create a series over the chunks
    data = collect(view(A, _chunk_inds(gc, I)...) for I in ci)
    return GeoSeries(data, DD.basedims(dims(A)))
end

# See iterate(::GridChunks) in Diskarrays.jl
function _chunk_inds(g, ichunk)
    outinds = map(ichunk.I, g.chunksize, g.parentsize, g.offset) do ic, cs, ps, of
        max((ic - 1) * cs + 1 -of, 1):min(ic * cs - of, ps)
    end
end

"""
    points(A::AbstractGeoArray; dims=(YDim, XDim))

Returns a generator of the points in `A` for dimensions in `dims`,
where points are a tuple of the values in each specified dimension
index.

The order of `dims` determines the order of the points.

$EXPERIMENTAL
"""
function points(A::AbstractGeoArray; dims=(YDim, XDim), ignore_missing=false)
    ignore_missing ? _points(A, dims) : _points_missing(A, dims)
end

function _points(A::AbstractGeoArray, dims)
    # Get the actual dims
    dims = DD.dims(A, dims)
    # Get the axes of each dimension
    dim_axes = map(d -> axes(d, 1), dims)
    # Construct a CartesianIndices generator
    indices = CartesianIndices(dim_axes)
    # Lazily index into the dimensions with the generator
    return (map(getindex, dims, Tuple(I)) for I in indices)
end
function _points_missing(A::AbstractGeoArray, dims)
    # Get the actual dims
    dims = DD.dims(A, dims)
    # Get the axes of each dimension
    dim_axes = map(d -> axes(d, 1), dims)
    # Construct a CartesianIndices generator
    indices = CartesianIndices(dim_axes)
    # Lazily index into the dimensions with the generator
    # or return missing if the matching array value is missing
    return (A[I] === missingval(A) ? map(getindex, dims, Tuple(I)) : missing for I in indices)
end

"""
	resample(A::AbstractGeoArray, resolution::Number; crs, method)
	resample(A::AbstractGeoArray; to::AbstractGeoArray, method)

`resample` uses `ArchGDAL.gdalwarp` to resample an `AbstractGeoArray`.

# Arguments

- `A`: The `AbstractGeoArray` to resample.
- `resolution`: A `Number` specifying the resolution for the output.
    If the keyword argument `crs` (described below) is specified, `resolution` must be in units of the `crs`.

# Keywords

- `to`: an `AbstractGeoArray` whos resolution, crs and bounds will be snapped to.
    For best results it should roughly cover the same extent, or a subset of `A`.
- `crs`: A `GeoFormatTypes.GeoFormat` specifying an output crs
    (`A` will be reprojected to `crs` in addition to being resampled). Defaults to `crs(A)`
- `method`: A `Symbol` or `String` specifying the method to use for resampling. Defaults to `:near`
    (nearest neighbor resampling). See [resampling method](https://gdal.org/programs/gdalwarp.html#cmdoption-gdalwarp-r)
    in the gdalwarp docs for a complete list of possible values.

$EXPERIMENTAL
"""
function resample end
function resample(A::AbstractGeoArray, resolution::Number;
    crs::GeoFormat=crs(A), method=:near
)
    wkt = convert(String, convert(WellKnownText, crs))
    flags = Dict(
        :t_srs => wkt,
        :tr => [resolution, resolution],
        :r => method,
    )
    return warp(A, flags)
end
function resample(A::AbstractGeoArray; to, method=:near)
    wkt = convert(String, convert(WellKnownText, crs(to)))
    latres, lonres = map(abs ∘ step, span(to, (Y(), X())))
    (latmin, latmax), (lonmin, lonmax) = bounds(to, (Y(), X()))
    flags = Dict(
        :t_srs => wkt,
        :tr => [latres, lonres],
        :te => [lonmin, latmin, lonmax, latmax],
        :r => method,
    )
    return warp(A, flags)
end
resample(st::AbstractGeoStack, args...; kw...) = map(A -> resample(A, args...; kw...), st)

"""
    warp(A::AbstractGeoArray, flags::Dict)

Gives access to the GDALs `gdalwarp` method given a `Dict` of flags,
where arguments than can be converted to strings, or vectors
of such arguments for flags that take multiple space separated arguments.

Arrays with additional dimensions not handled by GDAL (ie other than X, Y, Band)
are sliced, warped, and then combined - these dimensions will not change.

See: https://gdal.org/programs/gdalwarp.html for a list of arguments.

## Example

This simply resamples the array with the `:tr` (output file resolution) and `:r` flags:

```julia
using GeoData, RasterDataSources, Plots
A = GeoArray(WorldClim{Climate}, :prec; month=1)
flags = Dict(
    :tr => [1.0, 1.0],
    :r => :near,
)
warp(A, flags) |> plot
```

In practise, prefer [`resample`](@ref) for this. But `warp` may be more flexible.

$EXPERIMENTAL
"""
function warp(A::AbstractGeoArray, flags::Dict)
    odims = otherdims(A, (X, Y, Band))
    if length(odims) > 0
        # Handle dimensions other than X, Y, Band
        slices = slice(A, odims)
        warped = map(A -> _warp(A, flags), slices)
        return combine(warped, odims)
    else
        return _warp(A, flags)
    end
end
warp(st::AbstractGeoStack, flags::Dict) = map(A -> warp(A, flags), st)

function _warp(A::AbstractGeoArray, flags::Dict)
    flagvect = reduce([flags...]; init=[]) do acc, (key, val)
        append!(acc, String[_asflag(key), _stringvect(val)...])
    end
    AG.Dataset(A) do dataset
        AG.gdalwarp([dataset], flagvect) do warped
            _maybe_permute_from_gdal(read(GeoArray(warped)), dims(A))
        end
    end
end

_asflag(x) = string(x)[1] == '-' ? x : string("-", x)

_stringvect(x::AbstractVector) = Vector(string.(x))
_stringvect(x::Tuple) = [map(string, x)...]
_stringvect(x) = [string(x)]

"""
    mosaic(f, layers::AbstractGeoArray...; dims, missingval)
    mosaic(f, layers::AbstractGeoStack...; dims, missingval)
    mosaic(f, layers::Tuple; dims, missingval)

Combine `layer`s using the function `f`, (e.g. `mean`, `sum`,
`first` or `last`) where values from `layers` overlap.

# Keywords
- `dims`: The dimesions of `layers` to mosaic over, `(XDim, YDim)` by default.
    If dims contains an index it will be ignored, but this may change in future.
- `missingval`: Fills empty areas, and defualts to the
    `missingval` of the first layer.

$EXPERIMENTAL
"""
mosaic(f::Function, layers...; kw...) = mosaic(f, layers; kw...)
function mosaic(f::Function, layers::Tuple{<:AbstractGeoArray,Vararg};
    missingval=missingval(first(layers)), filename=nothing, kw...
)
    missingval isa Nothing && throw(ArgumentError("Layers have no missingval, so pass a `missingval` keyword explicitly"))
    T = Base.promote_type(typeof(missingval), Base.promote_eltype(layers...))
    dims = mosaic(map(DD.dims, layers))
    data = if filename isa Nothing
        Array{T,length(dims)}(undef, map(length, dims))
    else
        filename = create(filename, T, dims; missingval, metadata=metadata(first(layers)))
        parent(GeoArray(filename))
    end
    A = rebuild(first(layers); data, dims, missingval)
    open(A; write=true) do a
        mosaic!(f, a, layers; missingval, kw...)
    end
end
function mosaic(f::Function, layers::Tuple{<:AbstractGeoStack,Vararg}; kw...)
    map(layers...) do A...
        mosaic(f, A...; kw...)
    end
end

"""
    mosaic!(f, A, layers::AbstractGeoArray...; missingval, atol)
    mosaic!(f, A, layers::AbstractGeoStack...; missingval, atol)
    mosaic!(f, A, layers::Tuple; dims, missingval)

Combine `layer`s to the array `A` using the function `f`,
(e.g. `mean`, `sum`, `first` or `last`)
where values from `layers` overlap.

A may be a an opened disk-based `GeoArray` - the result will be
written to disk. `layers` should be memory-backed, or may experience
slow read speed with the current algorithm

# Keywords

- `missingval`: Fills empty areas, and defualts to the `missingval`
    of the first layer.
- `atol`: Absolute tolerance for comparison between index values.
    This is often required due to minor differences in range values
    due to floating point error. It is not applied to non-float dimensions.
    A tuple of tolerances may be passed, matching the dimension order.

$EXPERIMENTAL
"""
function mosaic!(f::Function, A, layers; missingval=missingval(A), atol=nothing)
    _without_mapped_crs(A) do A
        broadcast!(A, DimKeys(A; atol)) do ds
            # Get all the layers that have this point
            ls = foldl(layers; init=()) do acc, l
                if DD.hasselection(l, ds)
                    v = l[ds...]
                    (acc..., l)
                else
                    acc
                end
            end
            values = foldl(ls; init=()) do acc, l
                v = l[ds...]
                if isnothing(GeoData.missingval(l))
                    (acc..., v)
                elseif ismissing(GeoData.missingval(l))
                    ismissing(v) ? acc : (acc..., v)
                else
                    v === GeoData.missingval(l) ? acc : (acc..., v)
                end
            end
            if length(values) === 0
                missingval
            else
                f(values)
            end
        end
    end
end
function mosaic!(f::Function, A::AbstractGeoStack, stacks; kw...)
    map(A, stacks...) do a, s...
        mosaic!(f, a, s...; kw...)
    end
end

function mosaic(alldims::Tuple{<:DimTuple,Vararg{<:DimTuple}})
    map(mosaic, alldims...)
end
function mosaic(dims::Dimension...)
    map(dims) do d
        DD.comparedims(first(dims), d; val=false, length=false, mode=true)
    end
    _mosaic(mode(first(dims)), dims)
end

function _mosaic(mode::Categorical, dims::DimTuple)
    newindex = sort(union(map(val, dims)...); order=DD._ordering(indexorder(mode)))
    return rebuild(first(dims), newindex)
end
function _mosaic(mode::AbstractSampled, dims::DimTuple)
    order(mode) isa Unordered && throw(ArgumentError("Cant mozaic an Unordered dimension $(basetypeof(dim))"))
    _mosaic(span(mode), mode, dims)
end
function _mosaic(span::Regular, mode::AbstractSampled, dims::DimTuple)
    allkeys = map(val, dims)
    newindex = if indexorder(mode) isa ForwardIndex
        mi = minimum(map(first, allkeys))
        ma = maximum(map(last, allkeys))
        mi:step(span):ma
    else
        mi = minimum(map(last, allkeys))
        ma = maximum(map(first, allkeys))
        ma:step(span):mi
    end
    return rebuild(first(dims), newindex)
end
function _mosaic(::Irregular, mode::AbstractSampled, dims::DimTuple)
    newindex = sort(union(map(val, dims)...); order=DD._ordering(indexorder(mode)))
    return rebuild(first(dims), newindex)
end
function _mosaic(span::Explicit, mode::AbstractSampled, dims::DimTuple)
    newindex = sort(union(map(val, dims)...); order=DD._ordering(indexorder(mode)))
    bounds = map(val ∘ DD.span, dims)
    lower = map(b -> view(b, 1, :), bounds)
    upper = map(b -> view(b, 2, :), bounds)
    newlower = sort(union(lower...); order=DD._ordering(indexorder(mode)))
    newupper = sort(union(upper...); order=DD._ordering(indexorder(mode)))
    newbounds = vcat(permutedims(newlower), permutedims(newupper))
    newmode = rebuild(mode; span=Explicit(newbounds))
    return rebuild(first(dims); val=newindex, mode=newmode)
end

_without_mapped_crs(f, A) = _without_mapped_crs(f, A, mappedcrs(A))
_without_mapped_crs(f, A, ::Nothing) = f(A)
function _without_mapped_crs(f, A, mappedcrs)
    A = setmappedcrs(A, nothing)
    A = f(A)
    A = setmappedcrs(A, mappedcrs)
    return A
end

"""
   extract(vsr::AbstractGeoArray, points...; kw...)
   extract(vsr::AbstractGeoArray, points; kw...)

Extracts the value of the raster at the given points.
"""
function extract(A::AbstractGeoArray, x, y, others...; kw...)
    extract(A, (x, y, others...); kw...)
end
function extract(
    A::AbstractGeoArray, points; 
    order=(XDim, YDim, ZDim), atol=nothing
)
    ordereddims = dims(A, order)
    dimtypes = map(DD.basetypeof, ordereddims)
    seldims = map(ordereddims, dimtypes, points) do d, D, p
        if sampling(d) isa Points
            D(At(p; atol))
        else # Intervals
            D(Contains(p))
        end
    end
    if DD.hasselection(A, seldims)
        return A[seldims...]
    else
        return missing
    end
end
function extract(
    A::AbstractGeoArray, points::NTuple{N,T}; kw...
) where {T<:AbstractVector,N}
    extract.(Ref(A), zip(points...); kw...)
end
function extract(A::AbstractGeoArray, points::AbstractVector; kw...)
    extract.(Ref(A), points; kw...)
end
