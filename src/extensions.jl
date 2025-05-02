# extensions
function throw_extension_error(f::Function, package::String, extension::Symbol, args)
    @static if isdefined(Base, :get_extension) # julia > 1.9
    if isnothing(Base.get_extension(Rasters, extension))
        throw(BackendException(package))
    else
        throw(MethodError(f, args))
    end
    else
        throw(BackendException(package))
    end
end


# stubs that need ArchGDAL

"""
	resample(x; kw...)
    resample(xs...; to=first(xs), kw...)

`resample` uses `warp` (which uses GDALs `gdalwarp`) to resample a [`Raster`](@ref)
or [`RasterStack`](@ref) to a new `resolution` and optionally new `crs`,
or to snap to the bounds, resolution and crs of the object `to`.

Dimensions without an `AbstractProjected` lookup (such as a `Ti` dimension)
are iteratively resampled with GDAL and joined back into a single array.

If projections can be converted for each axis independently, it may
be faster and more accurate to use [`reproject`](@ref).

Run `using ArchGDAL` to make this method available.

# Arguments

- `x`: the object/s to resample.

# Keywords

- `to`: a `Raster`, `RasterStack`, `Tuple` of `Dimension` or `Extents.Extent`.
    If no `to` object is provided the extent will be calculated from `x`,
$RES_KEYWORD
$SIZE_KEYWORD
$CRS_KEYWORD
- `method`: A `Symbol` or `String` specifying the method to use for resampling.
    From the docs for [`gdalwarp`](https://gdal.org/programs/gdalwarp.html#cmdoption-gdalwarp-r):
    * `:near`: nearest neighbour resampling (default, fastest algorithm, worst interpolation quality).
    * `:bilinear`: bilinear resampling.
    * `:cubic`: cubic resampling.
    * `:cubicspline`: cubic spline resampling.
    * `:lanczos`: Lanczos windowed sinc resampling.
    * `:average`: average resampling, computes the weighted average of all non-NODATA contributing pixels.
        rms root mean square / quadratic mean of all non-NODATA contributing pixels (GDAL >= 3.3)
    * `:mode`: mode resampling, selects the value which appears most often of all the sampled points.
    * `:max`: maximum resampling, selects the maximum value from all non-NODATA contributing pixels.
    * `:min`: minimum resampling, selects the minimum value from all non-NODATA contributing pixels.
    * `:med`: median resampling, selects the median value of all non-NODATA contributing pixels.
    * `:q1`: first quartile resampling, selects the first quartile value of all non-NODATA contributing pixels.
    * `:q3`: third quartile resampling, selects the third quartile value of all non-NODATA contributing pixels.
    * `:sum`: compute the weighted sum of all non-NODATA contributing pixels (since GDAL 3.1)

    Where NODATA values are set to `missingval`.
$FILENAME_KEYWORD
$SUFFIX_KEYWORD

Note:
- GDAL may cause some unexpected changes in the raster, such as changing the `crs`
    type from `EPSG` to `WellKnownText` (it will represent the same CRS).

# Example

Resample a WorldClim layer to match an EarthEnv layer:

```jldoctest
using Rasters, RasterDataSources, ArchGDAL, Plots
A = Raster(WorldClim{Climate}, :prec; month=1)
B = Raster(EarthEnv{HabitatHeterogeneity}, :evenness)

a = plot(A)
b = plot(resample(A; to=B))

savefig(a, "build/resample_example_before.png");
savefig(b, "build/resample_example_after.png"); nothing

# output
```

### Before `resample`:

![before resample](resample_example_before.png)

### After `resample`:

![after resample](resample_example_after.png)

$EXPERIMENTAL
"""
resample(args...; kw...) = throw_extension_error(resample, "ArchGDAL", :RastersArchGDALExt, args)

"""
    warp(A::AbstractRaster, flags::Dict; kw...)

Gives access to the GDALs `gdalwarp` method given a `Dict` of
`flag => value` arguments that can be converted to strings, or vectors
where multiple space-separated arguments are required.

Arrays with additional dimensions not handled by GDAL (other than `X`, `Y`, `Band`)
are sliced, warped, and then combined to match the original array dimensions.
These slices will *not* be written to disk and loaded lazily at this stage -
you will need to do that manually if required.

See [the gdalwarp docs](https://gdal.org/programs/gdalwarp.html) for a list of arguments.

Run `using ArchGDAL` to make this method available.

# Keywords

$MISSINGVAL_KEYWORD
$FILENAME_KEYWORD
$SUFFIX_KEYWORD
- `missingval`: the missing value to use during warping, will default to
    `Rasters.missingval(A). Passing a pair will specify the missing value 
    to use after warping.

Any additional keywords are passed to `ArchGDAL.Dataset`.

## Example

This simply resamples the array with the `:tr` (output file resolution) and `:r`
flags, giving us a pixelated version:

```jldoctest
using Rasters, ArchGDAL, RasterDataSources, Plots
A = Raster(WorldClim{Climate}, :prec; month=1)
a = plot(A)

flags = Dict(
    :tr => [2.0, 2.0],
    :r => :near,
)
b = plot(warp(A, flags))

savefig(a, "build/warp_example_before.png");
savefig(b, "build/warp_example_after.png"); nothing

# output

```

### Before `warp`:

![before warp](warp_example_before.png)

### After `warp`:

![after warp](warp_example_after.png)

In practise, prefer [`resample`](@ref) for this. But `warp` may be more flexible.

$EXPERIMENTAL
"""
warp(args...; kw...) = throw_extension_error(warp, "ArchGDAL", :RastersArchGDALExt, args)

# stubs that need Proj

"""
    cellarea([method], x)

Gives the approximate area of each gridcell of `x`.
By assuming the earth is a sphere, it approximates the true size to about 0.1%, depending on latitude. 

Run `using ArchGDAL` or `using Proj` to make this method fully available.

- `method`: You can specify whether you want to compute the area in the plane of your projection `Planar()` or on a sphere of some radius `Spherical(; radius=...)`(the default).

- `Spherical` will compute cell area on the sphere, by transforming all points back to long-lat.  You can specify the radius by the `radius` keyword argument here.  By default, this is `6371008.8`, the mean radius of the Earth.

- `Planar` will compute cell area in the plane of the CRS you have chosen.  Be warned that this will likely be incorrect for non-equal-area projections.

Returns a Raster with the same x and y dimensions as the input, 
where each value in the raster encodes the area of the cell (in meters by default).

## Example

```julia
using Rasters, Proj, Rasters.Lookups
xdim = X(Projected(90.0:10.0:120; sampling=Intervals(Start()), crs=EPSG(4326)))
ydim = Y(Projected(0.0:10.0:50; sampling=Intervals(Start()), crs=EPSG(4326)))
myraster = rand(xdim, ydim)
cs = cellarea(myraster)

# output
╭───────────────────────╮
│ 4×6 Raster{Float64,2} │
├───────────────────────┴─────────────────────────────────────────────────── dims ┐
  ↓ X Projected{Float64} 90.0:10.0:120.0 ForwardOrdered Regular Intervals{Start},
  → Y Projected{Float64} 0.0:10.0:50.0 ForwardOrdered Regular Intervals{Start}
├───────────────────────────────────────────────────────────────────────── raster ┤
  extent: Extent(X = (90.0, 130.0), Y = (0.0, 60.0))

  crs: EPSG:4326
└─────────────────────────────────────────────────────────────────────────────────┘
   ↓ →  0.0        10.0        20.0        30.0            40.0      50.0
  90.0  1.23017e6   1.19279e6   1.11917e6   1.01154e6  873182.0  708290.0
 100.0  1.23017e6   1.19279e6   1.11917e6   1.01154e6  873182.0  708290.0
 110.0  1.23017e6   1.19279e6   1.11917e6   1.01154e6  873182.0  708290.0
 120.0  1.23017e6   1.19279e6   1.11917e6   1.01154e6  873182.0  708290.0
```
$EXPERIMENTAL
"""
cellarea(x; kw...) = cellarea(Spherical(), x; kw...)
cellarea(method::GeometryOpsCore.Manifold, x; kw...) = cellarea(method, dims(x, (XDim, YDim)); kw...)

function cellarea(method::GeometryOpsCore.Planar, dims::Tuple{<:XDim, <:YDim}; kw...)
    isintervals(dims) || throw(ArgumentError("Cannot calculate cell size for a `Raster` with `Points` sampling."))
    areas = _planar_cellarea(dims)
    return Raster(areas; dims)
end

function cellarea(method::GeometryOpsCore.Spherical, dims::Tuple{<:XDim, <:YDim}; kw...)
    isintervals(dims) || throw(ArgumentError("Cannot calculate cell size for a `Raster` with `Points` sampling."))
    areas = _spherical_cellarea(dims; radius = method.radius, kw...)
    return Raster(areas; dims)
end

_spherical_cellarea(args...; kw...) = throw_extension_error(_spherical_cellarea, "Proj", :RastersProjExt, args)

function _planar_cellarea(dims::Tuple{<:XDim, <:YDim})
    xbnds, ybnds = DD.intervalbounds(dims)
    broadcast(xb -> xb[2] - xb[1], xbnds) .* broadcast(yb -> yb[2] - yb[1], ybnds)'
end

function cellsize(args...; kw...) 
    @warn """
cellsize is deprecated and will be removed in a future version, use cellarea instead. 
Note that cellarea returns the area in square m, while cellsize still uses square km.
"""
    return cellarea(args...; kw..., radius = 6371.0088)
end

"""
Rasters.sample([rng], x, [n::Integer]; kw...)

Sample `n` random and optionally weighted points from from a `Raster` or `RasterStack`.
Returns a `Vector` of `NamedTuple`, closely resembling the return type of [`extract`](@ref).

Run `using StatsBase` to make this method available.
Note that this function is not exported to avoid confusion with StatsBase.sample

# Keywords

- `geometry`: include `:geometry` in returned `NamedTuple`. Specify the type and dimensions of the returned geometry by
    providing a `Tuple` or `NamedTuple` of dimensions. Defaults to `(X,Y)`
- `index`: include `:index` of the `CartesianIndex` in returned `NamedTuple`, `false` by default.
- `name`: a `Symbol` or `Tuple` of `Symbol` corresponding to layer/s of a `RasterStack` to extract. All layers by default.
- `skipmissing`: skip missing points automatically.
- `weights`: A DimArray that matches one or more of the dimensions of `x` with weights for sampling.
- `weightstype`: a `StatsBase.AbstractWeights` specifying the type of weights. Defaults to `StatsBase.Weights`.
- `replace`: sample with replacement, `true` by default. See `StatsBase.sample`
- `ordered`: sample in order, `false` by default. See `StatsBase.sample`

# Example
This code draws 5 random points from a raster, weighted by cell area.
```julia
using Rasters, Rasters.Lookups, Proj, StatsBase
xdim = X(Projected(90.0:10.0:120; sampling=Intervals(Start()), crs=EPSG(4326)))
ydim = Y(Projected(0.0:10.0:50; sampling=Intervals(Start()), crs=EPSG(4326)))
myraster = rand(xdim, ydim)
Rasters.sample(myraster, 5; weights=cellarea(myraster))

# output

5-element Vector{@NamedTuple{geometry::Tuple{Float64, Float64}, ::Union{Missing, Float64}}}:
 @NamedTuple{geometry::Tuple{Float64, Float64}, ::Union{Missing, Float64}}(((90.0, 10.0), 0.7360504790189618))
 @NamedTuple{geometry::Tuple{Float64, Float64}, ::Union{Missing, Float64}}(((90.0, 30.0), 0.5447657183842469))
 @NamedTuple{geometry::Tuple{Float64, Float64}, ::Union{Missing, Float64}}(((90.0, 30.0), 0.5447657183842469))
 @NamedTuple{geometry::Tuple{Float64, Float64}, ::Union{Missing, Float64}}(((90.0, 10.0), 0.7360504790189618))
 @NamedTuple{geometry::Tuple{Float64, Float64}, ::Union{Missing, Float64}}(((110.0, 10.0), 0.5291143028176258))
```
"""
sample(args...; kw...) = throw_extension_error(sample, "StatsBase", :RastersStatsBaseExt, args)



# Other shared stubs
function layerkeys end
function smapseries end
function maybe_correct_to_write end
function dims2geotransform end
function affine2geotransform end
function geotransform2affine end

# Shared between ArchGDAL and CoordinateTransformations extensions
const GDAL_EMPTY_TRANSFORM = [0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
const GDAL_TOPLEFT_X = 1
const GDAL_WE_RES = 2
const GDAL_ROT1 = 3
const GDAL_TOPLEFT_Y = 4
const GDAL_ROT2 = 5
const GDAL_NS_RES = 6

# The rest of the definition is in CoordinateTransformations
struct AffineProjected{T,F,A<:AbstractVector{T},M,C,MC,P,D} <: LA.Unaligned{T,1}
    affinemap::F
    data::A
    metadata::M
    crs::C
    mappedcrs::MC
    paired_lookup::P
    dim::D
end