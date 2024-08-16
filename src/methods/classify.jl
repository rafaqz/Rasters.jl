"""
    classify(x, pairs; lower=(>=), upper=(<), others=nothing)
    classify(x, pairs...; lower, upper, others)

Create a new array with values in `x` classified by the values in `pairs`.

`pairs` can hold tuples fo values `(2, 3)`, a `Fix2` function e.g. `<=(1)`, a `Tuple`
of `Fix2` e.g. `(>=(4), <(7))`, or an IntervalSets.jl interval, e.g. `3..9` or `OpenInterval(10, 12)`.
`pairs` can also be a `n * 3` matrix where each row is lower bounds, upper bounds, replacement.

If tuples or a `Matrix` are used, the `lower` and `upper` keywords define
how the lower and upper boundaries are chosen.

If `others` is set other values not covered in `pairs` will be set to that values.

# Arguments

- `x`: a `Raster` or `RasterStack`
- `pairs`: each pair contains a value and a replacement, a tuple of lower and upper
    range and a replacement, or a Tuple of `Fix2` like `(>(x), <(y)`.

# Keywords

- `lower`: Which comparison (`<` or `<=`) to use for lower values, if `Fix2` are not used.
- `upper`: Which comparison (`>` or `>=`) to use for upper values, if `Fix2` are not used.
- `others`: A value to assign to all values not included in `pairs`.
    Passing `nothing` (the default) will leave them unchanged.

# Example

```jldoctest
using Rasters, RasterDataSources, ArchGDAL, Plots
A = Raster(WorldClim{Climate}, :tavg; month=1)
classes = <=(15) => 10,
          15..25 => 20,
          25..35 => 30,
          >(35) => 40
classified = classify(A, classes; others=0, missingval=0)
plot(classified; c=:magma)

savefig("build/classify_example.png"); nothing

# output
```

![classify](classify_example.png)

$EXPERIMENTAL
"""
function classify end
classify(A::AbstractRaster, p1::Pair, pairs::Pair...; kw...) = classify(A, (p1, pairs...); kw...)
function classify(A::AbstractRaster, pairs::Union{Tuple,AbstractArray};
    filename=nothing, suffix=nothing, lower=(>=), upper=(<),
    others=nothing, missingval=missingval(A)
)
    # Make sure we get a concrete type. Broadcast doesn't always work.
    T = promote_type(_pairs_type(pairs), _others_type(others, A), typeof(missingval))
    # We use `Val{T}` to force type stability through the closure
    valT = Val{T}()
    f(x) = _convert_val(valT, _classify(x, pairs, lower, upper, others, Rasters.missingval(A), missingval))
    return create(filename, T, A; suffix, missingval) do C
        broadcast!(f, C, A)
    end
end
function classify(xs::AbstractRasterStack, pairs...; suffix=keys(xs), kw...)
    mapargs(xs, suffix) do x, s
        classify(x, pairs; suffix=s, kw...)
    end
end
function classify(xs::AbstractRasterSeries, pairs...; kw...)
    map(x -> classify(x, pairs...; kw...), xs)
end

_pairs_type(pairs) = promote_type(map(eltype ∘ last, pairs)...)
_pairs_type(pairs::AbstractArray{T}) where T = T

_others_type(others, A) = typeof(others)
_others_type(others::Nothing, A) = eltype(A)

_convert_val(::Val{T}, x) where T = convert(T, x)

"""
    classify!(x, pairs...; lower, upper, others)
    classify!(x, pairs; lower, upper, others)

Classify the values of `x` in-place, by the values in `pairs`.

If `Fix2` is not used, the `lower` and `upper` keywords

If `others` is set other values not covered in `pairs` will be set to that values.

# Arguments

- `x`: a `Raster` or `RasterStack`
- `pairs`: each pair contains a value and a replacement, a tuple of lower and upper
    range and a replacement, or a Tuple of `Fix2` like `(>(x), <(y)`.

# Keywords

- `lower`: Which comparison (`<` or `<=`) to use for lower values, if `Fix2` are not used.
- `upper`: Which comparison (`>` or `>=`) to use for upper values, if `Fix2` are not used.
- `others`: A value to assign to all values not included in `pairs`.
    Passing `nothing` (the default) will leave them unchanged.

# Example

`classify!` to disk, with key steps:
- copying a tempory file so we don't write over the RasterDataSources.jl version.
- use `open` with `write=true` to open the file with disk-write permissions.
- use `Float32` like `10.0f0` for all our replacement values and `other`, because
    the file is stored as `Float32`. Attempting to write some other type will fail.

```jldoctest
using Rasters, RasterDataSources, ArchGDAL, Plots
# Download and copy the file
filename = getraster(WorldClim{Climate}, :tavg; month=6)
tempfile = tempname() * ".tif"
cp(filename, tempfile)
# Define classes
classes = (5, 15) => 10,
          (15, 25) => 20,
          (25, 35) => 30,
          >=(35) => 40
# Open the file with write permission
open(Raster(tempfile); write=true) do A
    classify!(A, classes; others=0)
end
# Open it again to plot the changes
plot(Raster(tempfile); c=:magma)

savefig("build/classify_bang_example.png"); nothing

# output
```

![classify!](classify_bang_example.png)

$EXPERIMENTAL
"""
classify!(A::AbstractRaster, p1::Pair, pairs::Pair...; kw...) =
    classify!(A, (p1, pairs...); kw...)
function classify!(A::AbstractRaster, pairs;
    lower=(>=), upper=(<), others=nothing, missingval=missingval(A)
)
    T = promote_type(_pairs_type(pairs), _others_type(others, A), typeof(missingval))
    # We use `Val{T}` to force type stability through the closure
    valT = Val{T}()
    out = broadcast!(A, A) do x
        _convert_val(valT, _classify(x, pairs, lower, upper, others, Rasters.missingval(A), missingval))
    end
    return rebuild(out; missingval=missingval)
end
function classify!(xs::RasterSeriesOrStack, pairs...; kw...)
    map(x -> classify!(x, pairs...; kw...),  xs)
    return xs
end

# _classify
# Classify single values
function _classify(x, pairs, lower, upper, others, oldmissingval, newmissingval)
    isequal(x, oldmissingval) && return newmissingval
    # Use a fold instead of a loop, for type stability
    init = (false, last(first(pairs)))
    found, foundval = foldl(pairs; init) do (found, foundval), (find, replace)
        if !found && _compare(find, x, lower, upper)
            (true, replace)
        else
            (found, foundval)
        end
    end
    if found
        return foundval
    else
        return isnothing(others) ? x : others
    end
end
function _classify(x, pairs::AbstractArray, lower, upper, others, oldmissingval, newmissingval)
    isequal(x, oldmissingval) && return newmissingval
    found = false
    local foundval
    if ndims(pairs) != 2 || !(size(pairs, 2) in (2, 3))
        throw(ArgumentError("pairs must be a N*2 or N*3 matrix, or Pair"))
    elseif size(pairs, 2) == 2
        for i in 1:size(pairs, 1)
            find = pairs[i, 1]
            if _compare(find, x, lower, upper)
                foundval = pairs[i, 2]
                found = true
                break
            end
        end
    elseif size(pairs, 2) == 3
        for i in 1:size(pairs, 1)
            find = pairs[i, 1], pairs[i, 2]
            if _compare(find, x, lower, upper)
                foundval = pairs[i, 3]
                found = true
                break
            end
        end
    end
    if found
        return foundval
    else
        return isnothing(others) ? x : others
    end
end

_compare(find, x, lower, upper) = find === x
_compare(find::Base.Fix2, x, lower, upper) = find(x)
_compare((l, u)::Tuple, x, lower, upper) = lower(x, l) && upper(x, u)
_compare((l, u)::Tuple{<:Base.Fix2,<:Base.Fix2}, x, lower, upper) = l(x) && u(x)
_compare(interval::LA.IntervalSets.Interval, x, lower, upper) = x in interval
