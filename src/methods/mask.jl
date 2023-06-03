const TO_KEYWORD = """
- `to`: a `Raster`, `RasterStack`, `Tuple` of `Dimension` or `Extents.Extent`.
    If no `to` object is provided the extent will be calculated from the geometries,
    Additionally, when no `to` object or an `Extent` is passed for `to`, the `size`
    or `res` keyword must also be used.
"""
const SIZE_KEYWORD = """
- `size`: the size of the output array, as a `Tuple{Int,Int}` or single `Int` for a square.
    Only required when `to` is not used or is an `Extents.Extent`, and `res` is not used.
"""
const RES_KEYWORD = """
- `res`: the resolution of the dimensions, a `Real` or `Tuple{<:Real,<:Real}`.
    Only required when `to` is not used or is an `Extents.Extent`, and `size` is not used.
"""
const CRS_KEYWORD = """
- `crs`: a `crs` which will be attached to the resulting raster when `to` not passed
   or is an `Extent`. Otherwise the crs from `to` is used.
"""

const SHAPE_KEYWORDS = """
- `shape`: Force `data` to be treated as `:polygon`, `:line` or `:point` geometries.
    using points or lines as polygons may have unexpected results.
- `boundary`: for polygons, include pixels where the `:center` is inside the polygon,
    where the polygon `:touches` the pixel, or that are completely `:inside` the polygon.
    The default is `:center`.
"""

const GEOM_KEYWORDS = """
$TO_KEYWORD
$RES_KEYWORD
$SIZE_KEYWORD
$CRS_KEYWORD
$SHAPE_KEYWORDS
"""


"""
    mask(A:AbstractRaster; with, missingval=missingval(A))
    mask(x; with)

Return a new array with values of `A` masked by the missing values of `with`,
or by the shape of `with`, if `with` is a geometric object.

# Arguments

- `x`: a `Raster` or `RasterStack`

# Keywords

- `with`: an `AbstractRaster`, or any GeoInterface.jl compatible objects
    or table. The coordinate reference system of the point must match `crs(A)`.
- `missingval`: the missing value to use in the returned file.
- `filename`: a filename to write to directly, useful for large files.
- `suffix`: a string or value to append to the filename.
    A tuple of `suffix` will be applied to stack layers. `keys(st)` are the default.

# Geometry keywords

These can be used when `with` is a GeoInterface.jl compatible object:

$SHAPE_KEYWORDS

# Example

Mask an unmasked AWAP layer with a masked WorldClim layer,
by first resampling the mask.

```julia
using Rasters, RasterDataSources, ArchGDAL, Plots, Dates

# Load and plot the file
awap = read(Raster(AWAP, :tmax; date=DateTime(2001, 1, 1)))
a = plot(awap; clims=(10, 45))

# Create a mask my resampling a worldclim file
wc = Raster(WorldClim{Climate}, :prec; month=1)
wc_mask = resample(wc; to=awap)

# Mask
awap_masked = mask(awap; with=wc_mask)
b = plot(awap_masked; clims=(10, 45))

savefig(a, "build/mask_example_before.png");
savefig(b, "build/mask_example_after.png"); nothing
# output

```

### Before `mask`:

![before mask](mask_example_before.png)

### After `mask`:

![after mask](mask_example_after.png)

$EXPERIMENTAL
"""
mask(x; with, kw...) = _mask(x, with; kw...)

# Geometry mask
function _mask(s::AbstractRasterSeries, with; kw...)
    B = boolmask(with; to=dims(first(s), DEFAULT_POINT_ORDER), kw...)
    return _mask(s, B)
end
function _mask(x::RasterStackOrArray, with; kw...)
    # Geometries can only have `X`/`Y`/`Z` dims so limit them here
    B = boolmask(with; to=dims(x, DEFAULT_POINT_ORDER), kw...)
    return _mask(x, B)
end
# Array mask
function _mask(A::AbstractRaster, with::AbstractRaster;
    filename=nothing, suffix=nothing, missingval=_missingval_or_missing(A), kw...
)
    A1 = create(filename, A; suffix, missingval)
    open(A1; write=true) do a
        # The values array will be be written to A1 in `mask!`
        mask!(a; with, missingval, values=A)
    end
    return A1
end
function _mask(xs::AbstractRasterStack, with::AbstractRaster; suffix=keys(xs), kw...)
    mapargs((x, s) -> mask(x; with, suffix=s, kw...), xs, suffix)
end
function _mask(xs::AbstractRasterSeries, with::AbstractRaster; kw...)
    map(x -> mask(x; with, kw...), xs)
end

"""
    mask!(x; with, missingval=missingval(A))

Mask `A` by the missing values of `with`, or by all values outside `with` if it is a polygon.

If `with` is a polygon, creates a new array where points falling outside the polygon
have been replaced by `missingval(A)`.

Return a new array with values of `A` masked by the missing values of `with`,
or by a polygon.

# Arguments

- `x`: a `Raster` or `RasterStack`.

# Keywords

- `with`: another `AbstractRaster`, a `AbstractVector` of `Tuple` points,
    or any GeoInterface.jl `AbstractGeometry`. The coordinate reference system
    of the point must match `crs(A)`.
- `missingval`: the missing value to write to A in masked areas,
    by default `missingval(A)`.

# Example

Mask an unmasked AWAP layer with a masked WorldClim layer,
by first resampling the mask to match the size and projection.

```julia
using Rasters, RasterDataSources, ArchGDAL, Plots, Dates

# Load and plot the file
awap = read(RasterStack(AWAP, (:tmin, :tmax); date=DateTime(2001, 1, 1)))
a = plot(awap; clims=(10, 45), c=:imola)

# Create a mask my resampling a worldclim file
wc = Raster(WorldClim{Climate}, :prec; month=1)
wc_mask = resample(wc; to=awap)

# Mask
mask!(awap; with=wc_mask)
b = plot(awap; clims=(10, 45))

savefig(a, "build/mask_bang_example_before.png");
savefig(b, "build/mask_bang_example_after.png"); nothing

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
mask!(xs::RasterStackOrArray; with, kw...) = _mask!(xs, with; kw...)

# Geometry mask
function _mask!(x::RasterStackOrArray, geom; kw...)
    B = boolmask(geom; to=dims(x), kw...)
    _mask!(x, B; missingval=false)
    return x
end
# Array mask
function _mask!(st::RasterStack, with::AbstractRaster; kw...)
    map(A -> mask!(A; with, kw...), st)
    return st
end
function _mask!(A::AbstractRaster, with::AbstractRaster;
    missingval=missingval(A), values=A
)
    missingval isa Nothing && _nomissingerror()
    missingval = convert(eltype(A), missingval)

    broadcast_dims!(A, values, with) do x, w
        if isequal(w, Rasters.missingval(with))
            missingval
        else
            convert(eltype(A), x)
        end
    end
    return A
end

_nomissingerror() = throw(ArgumentError("Array has no `missingval`. Pass a `missingval` keyword compatible with the type, or use `rebuild(A; missingval=somemissingval)` to set it."))

"""
    boolmask(obj::Raster; [missingval])
    boolmask(obj; [to, res, size])

Create a mask array of `Bool` values, from another `Raster`. An
`AbstractRasterStack` or `AbstractRasterSeries` are also accepted, but a mask
is taken of the first layer or object *not* all of them.

The array returned from calling `boolmask` on a `AbstractRaster` is a
[`Raster`](@ref) with the same dimensions as the original array and a
`missingval` of `false`.

# Arguments

- `obj`: a [`Raster`](@ref), a GeoInterface.jl geometry, or a vector or table of geometries.

# `Raster` / `RasterStack` Keywords

- `missingval`: The missing value of the source array, with default `missingval(raster)`.

# Geometry keywords

$GEOM_KEYWORDS

And specifically for `shape=:polygon`:

- `boundary`: include pixels where the `:center` is inside the polygon, where
    the line `:touches` the pixel, or that are completely `:inside` inside the polygon.
    The default is `:center`.

For tabular data, feature collections and other iterables

- `collapse`: if `true`, collapse all geometry masks into a single mask. Otherwise
    return a Raster with an additional `geometry` dimension, so that each slice
    along this axis is the mask of the `geometry` opbject of each row of the
    table, feature in the feature collection, or just each geometry in the iterable.

# Example

```jldoctest
using Rasters, RasterDataSources, ArchGDAL, Plots, Dates
wc = Raster(WorldClim{Climate}, :prec; month=1)
boolmask(wc) |> plot

savefig("build/boolmask_example.png"); nothing

# output
```

![boolmask](boolmask_example.png)

$EXPERIMENTAL
"""
function boolmask end
boolmask(series::AbstractRasterSeries; kw...) = boolmask(first(series); kw...)
boolmask(stack::AbstractRasterStack; kw...) = boolmask(first(stack); kw...)
function boolmask(source::AbstractRaster; kw...)
    dest = _init_bools(source, Bool, nothing; kw..., missingval=false)
    return boolmask!(dest, source; kw...)
end
function boolmask(x; to=nothing, kw...)
    if to isa Union{AbstractDimArray,AbstractDimStack,DimTuple}
        to = dims(to, DEFAULT_POINT_ORDER)
    end
    A = _init_bools(to, Bool, x; kw..., missingval=false)
    return boolmask!(A, x; kw...)
end

function boolmask!(dest::AbstractRaster, src::AbstractRaster;
    missingval=_missingval_or_missing(src)
)
    broadcast_dims!(dest, src) do a
        !isequal(a, missingval)
    end
end
function boolmask!(dest::AbstractRaster, geoms;
    allocs=nothing, lock=nothing, progress=true, threaded=true, kw...
)
    if isnothing(allocs)
        allocs = _burning_allocs(dest; threaded)
    end
    if hasdim(dest, :geometry)
        range = _geomindices(geoms)
        _run(range, threaded, progress, "Burning each geometry to a BitArray slice...") do i
            geom = _getgeom(geoms, i)
            ismissing(geom) && return nothing
            slice = view(dest, Dim{:geometry}(i))
            # We don't need locks - these are independent slices
            burn_geometry!(slice, geom; kw..., fill=true, allocs=_get_alloc(allocs))
            return nothing
        end
    else
        burn_geometry!(dest, geoms; kw..., allocs, lock, progress, threaded, fill=true)
    end
    return dest
end

"""
    missingmask(obj::Raster; kw...)
    missingmask(obj; [to, res, size, collapse])

Create a mask array of `missing` and `true` values, from another `Raster`.
`AbstractRasterStack` or `AbstractRasterSeries` are also accepted, but a mask
is taken of the first layer or object *not* all of them.

For [`AbstractRaster`](@ref) the default `missingval` is `missingval(A)`,
but others can be chosen manually.

The array returned from calling `missingmask` on a `AbstractRaster` is a
[`Raster`](@ref) with the same size and fields as the original array.

# Keywords

$GEOM_KEYWORDS

# Example

```jldoctest
using Rasters, RasterDataSources, ArchGDAL, Plots, Dates
wc = Raster(WorldClim{Climate}, :prec; month=1)
missingmask(wc) |> plot

savefig("build/missingmask_example.png"); nothing

# output
```

![missingmask](missingmask_example.png)

$EXPERIMENTAL
"""
missingmask(series::AbstractRasterSeries; kw...) = missingmask(first(series); kw...)
missingmask(stack::AbstractRasterStack; kw...) = missingmask(first(stack); kw...)
function missingmask(source::AbstractRaster; kw...)
    dest = _init_bools(source, Union{Missing,Bool}, nothing; kw..., missingval=missing)
    return missingmask!(dest, source; kw...)
end
function missingmask(x; to=nothing, kw...)
    B = _init_bools(to, Union{Missing,Bool}, x; kw..., missingval=missing)
    return missingmask!(B, x; kw...)
end

function missingmask!(dest::AbstractRaster, src::AbstractRaster;
    missingval=_missingval_or_missing(src)
)
    broadcast_dims!(dest, src) do x
        isequal(x, missingval) ? missing : true
    end
end
function missingmask!(dest::AbstractRaster, geom; kw...)
    B = boolmask!(dest, geom; kw...)
    dest .= (b -> b ? true : missing).(B)
    return dest
end
