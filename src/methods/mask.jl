"""
    mask(A:AbstractRaster; to, missingval=missingval(A))
    mask(x; to, order=(XDim, YDim))

Return a new array with values of `A` masked by the missing values of `to`,
or by when more than 50% outside `to`, if it is a polygon.

# Arguments

- `x`: a `Raster` or `RasterStack`

# Keywords

- `to`: another `AbstractRaster`, a `AbstractVector` of `Tuple` points,
    or any GeoInterface.jl `AbstractGeometry`. The coordinate reference system
    of the point must match `crs(A)`.
- `order`: the order of `Dimension`s in the points. Defaults to `(XDim, YDim)`.
- `missingval`: the missing value to use in the returned file.
- `filename`: a filename to write to directly, useful for large files.
- `suffix`: a string or value to append to the filename.
    A tuple of `suffix` will be applied to stack layers. `keys(st)` are the default.

# Geometry keywords

These can be used when a `GeoInterface.AbstractGeometry` is passed in.

- `shape`: Force `data` to be treated as `:polygon`, `:line` or `:point`.

And specifically for `shape=:polygon`:

- `boundary`: include pixels where the `:center` is inside the polygon, where 
    the line `:touches` the pixel, or that are completely `:inside` inside the polygon.

In future this method will accept more point types.

# Example

Mask an unmasked AWAP layer with a masked WorldClim layer,
by first resampling the mask.

```jldoctest
using Rasters, Plots, Dates

# Load and plot the file
awap = read(Raster(AWAP, :tmax; date=DateTime(2001, 1, 1)))
a = plot(awap; clims=(10, 45))

# Create a mask my resampling a worldclim file
wc = Raster(WorldClim{Climate}, :prec; month=1)
wc_mask = resample(wc; to=awap)

# Mask
awap_masked = mask(awap; to=wc_mask)
b = plot(awap_masked; clims=(10, 45))

savefig(a, "build/mask_example_before.png")
savefig(b, "build/mask_example_after.png")
# output

```

### Before `mask`:

![before mask](mask_example_before.png)

### After `mask`:

![after mask](mask_example_after.png)

$EXPERIMENTAL
"""
mask(x; to, kw...) = _mask(x, to; kw...)

# Geometry mask
function _mask(s::AbstractRasterSeries, geom; kw...)
    B = boolmask(geom; template=first(s), kw...)
    return _mask(s, B)
end
function _mask(x::RasterStackOrArray, geom; kw...)
    B = boolmask(geom; template=x, kw...)
    return _mask(x, B)
end
# Array mask
function _mask(A::AbstractRaster, to::AbstractRaster; 
    filename=nothing, suffix=nothing, missingval=_missingval_or_missing(A), kw...
)
    A1 = create(filename, A; suffix, missingval)
    open(A1; write=true) do a
        # The values array will be be written to A1 in `mask!`
        mask!(a; to, missingval, values=A)
    end
    return A1
end
function _mask(xs::AbstractRasterStack, to::AbstractRaster; suffix=keys(xs), kw...) 
    mapargs((x, s) -> mask(x; to, suffix=s, kw...), xs, suffix)
end
function _mask(xs::AbstractRasterSeries, to::AbstractRaster; kw...) 
    map(x -> mask(x; to, kw...), xs)
end

"""
    mask!(A; to, missingval=missingval(A), order=(XDim, YDim))

Mask `A` by the missing values of `to`, or by values outside `to` if i is a polygon.

If `to` is a polygon, creates a new array where points falling outside the polygon
have been replaced by `missingval(A)`.

Return a new array with values of `A` masked by the missing values of `to`,
or by a polygon.

# Arguments

- `x`: a `Raster` or `RasterStack`.

# Keywords

- `to`: another `AbstractRaster`, a `AbstractVector` of `Tuple` points,
    or any GeoInterface.jl `AbstractGeometry`. The coordinate reference system
    of the point must match `crs(A)`.
- `order`: the order of `Dimension`s in the points. Defaults to `(XDim, YDim)`.
- `missingval`: the missing value to write to A in masked areas,
    by default `missingval(A)`.

# Example

Mask an unmasked AWAP layer with a masked WorldClim layer,
by first resampling the mask to match the size and projection.

```jldoctest
using Rasters, Plots, Dates

# Load and plot the file
awap = read(RasterStack(AWAP, (:tmin, :tmax); date=DateTime(2001, 1, 1)))
a = plot(awap; clims=(10, 45))

# Create a mask my resampling a worldclim file wc = Raster(WorldClim{Climate}, :prec; month=1)
wc_mask = resample(wc; to=awap)
# Mask mask!(awap; to=wc_mask) b = plot(awap; clims=(10, 45))

savefig(a, "build/mask_bang_example_before.png")
savefig(b, "build/mask_bang_example_after.png")
# output
```

### Before `mask!`:

![before mask!](mask_bang_example_before.png)

### After `mask!`:

![after mask!](mask_bang_example_after.png)

$EXPERIMENTAL
"""
function mask! end
function mask!(xs::AbstractRasterSeries; kw...)
    foreach(x -> mask!(x; kw...), xs)
    return xs
end
mask!(xs::RasterStackOrArray; to, kw...) = _mask!(xs, to; kw...)

# Geometry mask
function _mask!(x::RasterStackOrArray, geom; kw...)
    B = boolmask(geom; template=dims(x), kw...)
    _mask!(x, B; missingval=fals)
    return x
end
# Array mask
function _mask!(st::RasterStack, to::AbstractRaster; kw...)
    map(A -> mask!(A; to, kw...), st)
    return st
end
function _mask!(A::AbstractRaster, to::AbstractRaster;
    missingval=missingval(A), values=A
)
    missingval isa Nothing && _nomissingerror()
    missingval = convert(eltype(A), missingval)

    broadcast_dims!(A, values, to) do s, t
        isequal(t, Rasters.missingval(to)) ? missingval : convert(eltype(A), s)
    end
    return A
end

_nomissingerror() = throw(ArgumentError("Array has no `missingval`. Pass a `missingval` keyword compatible with the type, or use `rebuild(A; missingval=somemissingval)` to set it."))

"""
    boolmask(A::AbstractArray; [missingval])
    boolmask(geom; template, order, shape, crossing)

Create a mask array of `Bool` values, from any `AbstractArray`.

The array returned from calling `boolmask` on a `AbstractRaster` is a
[`Raster`](@ref) with the same size and fields as the original array.

# Arguments

- `T`: `BitArray` or `Array`
- `A`: An `AbstractArray`.

# Keywords

For arrays:

- `missingval`: The missing value of the source array. For [`AbstractRaster`](@ref) the
    default `missingval` is `missingval(A)`, for all other `AbstractArray`s it is `missing`.

For gemoetries:

- `template`: an `AbstractRaster` or `AbstractRasterStack`.
- `order`: the order of `Dimension`s in the points. Defaults to `(XDim, YDim)`.

# Example

```jldoctest
using Rasters, Plots, Dates
wc = Raster(WorldClim{Climate}, :prec; month=1)
boolmask(wc) |> plot

savefig("build/boolmask_example.png")
# output
```

![boolmask](boolmask_example.png)

$EXPERIMENTAL
"""
function boolmask end
boolmask(x::RasterStackOrArray; kw...) = boolmask!(_bools(x, dims(x)), x; kw...)
function boolmask(x; template, order=(XDim, YDim), kw...)
    boolmask!(_bools(template, commondims(template, order)), x; order, kw...)
end

_bools(x) = _bools(x, dims(x))
_bools(x::AbstractRasterSeries, dims) = _bools(first(x), dims)
_bools(x::AbstractRasterStack, dims) = _bools(first(x), dims)
_bools(x, dims) = Raster(falses(dims); missingval=false)
function _bools(x::AbstractRaster, dims)
    data = if parent(x) isa Array || parent(x) isa DA.AbstractDiskArray
        falses(dims) # Use a BitArray
    else
        fill!(similar(parent(x), Bool), false)
    end
    return Raster(data; dims, missingval=false)
end

function boolmask!(dest::AbstractRaster, src::AbstractRaster; 
    missingval=_missingval_or_missing(src)
)
    broadcast!(a -> !isequal(a, missingval), dest, src)
end
boolmask!(dest::AbstractRaster, geom; kw...) = _fill_geometry!(dest, geom; kw...)

"""
    missingmask(x; kw...)

Create a mask array of `missing` or `true` values, from any `AbstractArray`.
For [`AbstractRaster`](@ref) the default `missingval` is `missingval(A)`,
for all other `AbstractArray`s it is `missing`.

The array returned from calling `missingmask` on a `AbstractRaster` is a
[`Raster`](@ref) with the same size and fields as the original array.

# Example

```jldoctest
using Rasters, Plots, Dates
wc = Raster(WorldClim{Climate}, :prec; month=1)
missingmask(wc) |> plot

savefig("build/missingmask_example.png")
# output
```

![missingmask](missingmask_example.png)

$EXPERIMENTAL
"""
function missingmask(x::AbstractRaster; missingval=missingval(x))
    B = Raster(zeros(Union{Bool,Missing}, dims(x)); missingval=missing)
    B .= missing
    boolmask!(B, x; missingval)
    return broadcast(x -> x ? true : missing, B)
end
