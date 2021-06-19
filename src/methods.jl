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
"""
function crop end
crop(layers::AbstractGeoArray...; kw...) = crop(layers)
function crop(layers::Union{Tuple,NamedTuple}; to=_smallestdims(layers))
    map(l -> crop(l; to), layers)
end
crop(A::AbstractGeoArray; to) = _cropto(A, to)

# crop `A` to values of dims of `to`
_cropto(A::AbstractGeoArray, to) = _cropto(A, dims(to)) 
function _cropto(A::AbstractGeoArray, to::Tuple)
    # Create selectors for each dimension
    # `Between` the bounds of the dimension
    selectors = map(to) do d
        DD.basetypeof(d)(Between(DD.bounds(d)))
    end
    # Take a view of the `Between` selectors
    return view(A, selectors...)
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
    extend(A::AbstractGeoArray; to::Tuple)

Extend multiple [`AbstractGeoArray`](@ref) to match
the size of the largest one. A single `AbstractGeoArray` 
can be extended bu passing the new `dims` tuple as the second
argumen.
"""
function extend end
extend(layers::AbstractGeoArray...) = extend(layers)
function extend(layers::Union{NamedTuple,Tuple}; to=_largestdims(layers))
    # Extend all layers to `to`, by default the _largestdims
    map(l -> extend(l; to), layers)
end
function extend(A::AbstractGeoArray; to::Tuple)
    size = map(length, to)
    T = eltype(A)
    # Creat a new extended array
    newdata = similar(parent(A), T, size)
    # Fill it with missing/nodata values
    newdata .= missingval(A)
    # Rebuild the original object with larger data and dims.
    newA = rebuild(A; data=newdata, dims=to)
    # Calculate the range of the old array in the extended array
    ranges = map(dims(A), to) do d, nd
        (DD.sel2indices(nd, Near(first(d)))):(DD.sel2indices(nd, Near(last(d))))
    end
    # Copy the original data to the new array
    copyto!(
        parent(newA), CartesianIndices((ranges...,)), 
        parent(read(A)), CartesianIndices(A)
    ) 
    return newA
end

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

The trimmed size will be padded by `pad` on all sides.
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

# Get the ranges to trim to for dimensions in `dims`
function _trimranges(A, targetdims)
    # Bool vectors to track wich rows/cols have non-missing values
    # for the dimensions we are interested in
    trackers = map(s -> zeros(Bool, s), map(d -> size(A, d), targetdims))
    # Broadcast over the array and index generators
    index_generators = DD.dimwise_generators(dims(A))
    updates = broadcast(A, index_generators) do a, I
        # Check if the value is non-missing
        if a !== missingval(A)
            # Set the tracker for this index to true. We are only tracking
            # the target dims, so `dims` extracts them from the tuple I.
            inds = map(val, dims(I, targetdims))
            for (i, n) in enumerate(inds)
                trackers[i][n] = true
            end
        end
        nothing
    end
    collect(updates)
    # Get the ranges that contain all non-missing values
    cropranges = map(trackers) do t
        findfirst(t):findlast(t)
    end
    return cropranges
end

"""
    slice(A::Union{AbstractGeoArray,AbstractGeoStack,AbstracGeoSeries}, dims)

Slice an object allong some dimension(s), lazily using `view`. For a single `GeoArray` 
or `GeoStack` this will return a `GeoSeries` of `GeoArray` or `GeoStack` that are slices 
along the specified dimensions. For a `GeoSeries`, the output is another series where
the child objects are sliced and the series dimensions index is now of the child 
dimensions combined. `slice` on a `GeoSeries` with no dimensions will slice along
the dimensions shared by both the series and child object.
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

