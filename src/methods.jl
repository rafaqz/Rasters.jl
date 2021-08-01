"""
    replace_missing(a::AbstractGeoArray, newmissingval)
    replace_missing(a::AbstractGeoStack, newmissingval)

Replace missing values in the array or stack with a new missing value,
also updating the `missingval` field/s.

A `GeoArray` containing a newly allocated `Array` is always returned,
even when the missing value matches the current value.
"""
function replace_missing(A::AbstractGeoArray, newmissingval=missing)
    A = read(A)
    newdata = if ismissing(missingval(A))
        if ismissing(newmissingval)
            copy(parent(read(A)))
        else
            collect(Missings.replace(parent(A), newmissingval))
        end
    else
        replace(parent(A), missingval(A) => newmissingval)
    end
    rebuild(A; data=newdata, missingval=newmissingval)
end
function replace_missing(stack::AbstractGeoStack, newmissingval=missing)
    rebuild(stack, map(a -> replace_missing(a, newmissingval, values(stack))))
end

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
    crop(layers::AbstractGeoArray...)
    crop(layers::Union{NamedTuple,Tuple})
    crop(A::AbstractGeoArray; to::Tuple)

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
crop(A::AbstractGeoArray; to) = _crop_to(A, to)

# crop `A` to values of dims of `to`
_crop_to(A::AbstractGeoArray, to) = _crop_to(A, dims(to)) 
function _crop_to(A::AbstractGeoArray, to::Tuple)
    # Create selectors for each dimension
    # `Between` the bounds of the dimension
    selectors = map(to) do d
        DD.basetypeof(d)(Between(DD.bounds(d)))
    end
    # Take a view of the `Between` selectors
    return _without_mapped_crs(A) do a
        view(a, selectors...)
    end
end

_without_mapped_crs(f, A) = _without_mapped_crs(f, A, mappedcrs(A))
_without_mapped_crs(f, A, ::Nothing) = f(A)
function _without_mapped_crs(f, A, mappedcrs)
    # Drop mappedcrs
    A = set(A; 
        X=rebuild(mode(A, X); mappedcrs=nothing),
        Y=rebuild(mode(A, Y); mappedcrs=nothing),
    )
    A = f(A)
    # Re-apply mappedcrs 
    A = set(A; 
        X=rebuild(mode(A, X); mappedcrs=mappedcrs),
        Y=rebuild(mode(A, Y); mappedcrs=mappedcrs),
    )
    return A
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
    extend(A::AbstractGeoArray; to)

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
extend(A::AbstractGeoArray; to=dims(A)) = _extend_to(A, to)

_extend_to(A::AbstractGeoArray, to) = _extend_to(A, dims(to))
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

@inline _maybeflipindex(::ForwardRelation, d, i) = i
@inline _maybeflipindex(::ReverseRelation, d, i::Integer) = lastindex(d) - i + 1
@inline _maybeflipindex(::ReverseRelation, d, i::AbstractArray) =
    reverse(lastindex(d) .- i .+ 1)

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
        
Trim `missingval` from `A` for axes in dims.

By default `dims=(X, Y)`, so trimming keeps the area of `X` and `Y` 
that contains non-missing values along all other dimensions.

The trimmed size will be padded by `pad` on all sides, although 
padding will not be added beyond the original extent of the array.

$EXPERIMENTAL
"""
function trim(A::GeoArray; dims::Tuple=(X(), Y()), pad::Int=0)
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
    # Wrap indices in diensions so we can sort and filter them
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
    trackers .= A .!== missingval(A)
    # Get the ranges that contain all non-missing values
    cropranges = map(trackers.tracking) do a
        findfirst(a):findlast(a)
    end
    return cropranges
end

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
slice(x::Union{AbstractGeoArray,AbstractGeoStack}, dims) = slice(x, (dims,))
function slice(x::Union{AbstractGeoArray,AbstractGeoStack}, dims::Tuple)
    # Make sure all dimensions in `dims` are in `x`
    all(hasdim(x, dims)) || _errordimsnotfound(otherdims(dims, DD.dims(x)))
    # Define dimensions and data for the sliced GeoSeries
    seriesdims = DD.dims(x, dims)
    # series data is a generator of view slices
    seriesdata = map(DD.dimwise_generators(seriesdims)) do ds
        view(x, ds...)
    end
    return GeoSeries(seriesdata, seriesdims)
end
slice(ser::AbstractGeoSeries, dims) = cat(map(x -> slice(x, dims), ser)...; dims=dims)

@noinline _errordimsnotfound(dims) = 
    throw(ArgumentError("$(map(DD.dim2key, dims)) were dims not found in $(nameof(typeof(x)))"))
@noinline _errordimsnotfound(dims::Tuple{<:Any}) = 
    throw(ArgumentError("$(DD.dim2key(dims[1])) was dim not found in $(nameof(typeof(x)))"))

# Get the bounds wrapped in Dim(Between)
dimbounds(A::AbstractDimArray) = dimbounds(bounds, A)
function dimbounds(f::Function, A::AbstractDimArray)
    map(dims(A), f(A)) do dim, bounds
        # TODO is Between the right selector? do we need inclusive?
        basetypeof(dim)(Between(bounds))
    end
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
    points(A::AbstractGeoArray; dims=(Y, X))
    
Returns a generator of the points in `A` for dimensions in `dims`,
where points are a tuple of the values in each specified dimension 
index.

The order of `dims` determines the order of the points.

$EXPERIMENTAL
"""
function points(A::AbstractGeoArray; dims=(Y, X))
    # Get the actual dims
    dims = DD.dims(A, dims)
    # Get the axes of each dimension
    dim_axes = map(d -> axes(d, 1), dims)
    # Construct a CartesianIndices generator
    indices = CartesianIndices(dim_axes)
    # Lazily index into the dimensions with the generator
    return (map(getindex, dims, Tuple(I)) for I in indices)
end

