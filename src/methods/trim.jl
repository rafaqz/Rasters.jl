"""
    trim(x; dims::Tuple, pad::Int)

Trim `missingval(x)` from `x` for axes in `dims`, returning a view of `x`.

# Arguments

- `x`: A `Raster` or `RasterStack`. For stacks, all layers must having
    missing values for a pixel for it to be trimmed.

# Keywords

- `dims`: By default `dims=(XDim, YDim)`, so that trimming keeps the area 
    of `X` and `Y` that contains non-missing values along all other dimensions.
- `pad`: The trimmed size will be padded by `pad` on all sides, although
    padding will not be added beyond the original extent of the array.

`trim` does not accept `filename`/`suffix` arguments as it does not alter the underlying data.

# Example

Create trimmed layers of Australian habitat heterogeneity.

```jldoctest
using Rasters, RasterDataSources, Plots
layers = (:evenness, :range, :contrast, :correlation)
st = RasterStack(EarthEnv{HabitatHeterogeneity}, layers)

# Roughly cut out australia
ausbounds = X(100 .. 160), Y(-50 .. -10)
aus = st[ausbounds...]
a = plot(aus)

# Trim missing values and plot
b = plot(trim(aus))

savefig(a, "build/trim_example_before.png");
savefig(b, "build/trim_example_after.png"); nothing

# output

```

### Before `trim`:

![before trim](trim_example_before.png)

### After `trim`:

![after trim](trim_example_after.png)

$EXPERIMENTAL
"""
function trim(x::RasterStackOrArray; dims::Tuple=(XDim, YDim), pad::Int=0)
    # Get the actual dimensions in their order in the array
    targetdims = commondims(x, dims)
    # Get the range of non-missing values for each dimension
    trackers = AxisTrackers(DD.dims(x), targetdims)
    _update!(trackers, x)
    ranges = _cropranges(trackers)
    padded_ranges = _pad(ranges, targetdims, pad)
    rangedims = map(rebuild, targetdims, padded_ranges)
    return view(x, rangedims...)
end
function trim(ser::AbstractRasterSeries; dims::Tuple=(XDim, YDim), pad::Int=0)
    x = first(ser)
    # Get the actual dimensions in their order in the array
    targetdims = commondims(x, dims)
    # Get the range of non-missing values for each dimension
    trackers = AxisTrackers(DD.dims(x), targetdims)
    map(x -> _update!(trackers, x), ser)
    # Get the ranges that contain all non-missing values
    ranges = _cropranges(trackers)
    padded_ranges = _pad(ranges, targetdims, pad)
    rangedims = map(rebuild, targetdims, padded_ranges)
    return map(x -> view(x, rangedims...), ser)
end

function _pad(ranges, dims, pad)
    map(ranges, size(dims)) do r, l
        max(first(r)-pad, 1):min(last(r)+pad, l)
    end
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

function _cropranges(trackers::AxisTrackers)
    map(trackers.tracking) do a
        f = findfirst(a)
        l = findlast(a)
        f = f === nothing ? firstindex(a) : f
        l = l === nothing ? lastindex(a) : l
        f:l
    end
end

# Broadcast over the array and tracker to mark axis indices as being missing or not
_update!(tr::AxisTrackers, A::AbstractRaster) = tr .= A .!== missingval(A)
_update!(tr::AxisTrackers, st::AbstractRasterStack) = map(A -> tr .= A .!== missingval(A), st)

