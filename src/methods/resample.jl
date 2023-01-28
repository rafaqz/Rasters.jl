"""
	resample(x, resolution::Number; crs, method)
	resample(x; to, method)
    resample(xs...; to=first(xs), method)

`resample` uses `ArchGDAL.gdalwarp` to resample a [`Raster`](@ref) or
[`RasterStack`](@ref) to a new `resolution` and optionally new `crs`,
or to snap to the bounds, resolution and crs of the object `to`.

# Arguments

- `x`: the object to resample.
- `resolution`: a `Number` specifying the resolution for the output.
    If the keyword argument `crs` (described below) is specified, `resolution` must be in units of the `crs`.

# Keywords

- `to`: an `AbstractRaster` whose resolution, crs and bounds will be snapped to.
    For best results it should roughly cover the same extent, or a subset of `A`.
- `crs`: A `GeoFormatTypes.GeoFormat` such as `EPSG(x)` or `WellKnownText(string)` specifying an
    output crs (`A` will be reprojected to `crs` in addition to being resampled). Defaults to `crs(A)`
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

Note: `missingval` of `missing` does not work with GDAL. Use `replace_missing(A, newmissingval)` to
assign a missing value before using `resample` if the current value is `missing`. This will be
automated in future versions.

# Example

Resample a WorldClim layer to match an EarthEnv layer:

```jldoctest
using Rasters
using CairoMakie
CairoMakie.activate!()

A = Raster(WorldClim{Climate}, :prec; month=1)
B = Raster(EarthEnv{HabitatHeterogeneity}, :evenness)

#a = plot(A)
a = lines(1:10)
#b = plot(resample(A; to=B))
b = lines(1:10)
save("docs/assets/resample_example_before.png", a); # hide
save("docs/assets/resample_example_after.png", b); nothing
# output
```

### Before `resample`:

![before resample](./assets/resample_example_before.png) # hide

### After `resample`:

![after resample](./assets/resample_example_after.png) # hide

$EXPERIMENTAL
"""
function resample end
resample(xs::RasterStackOrArray...; kw...) = resample(xs; kw...)
function resample(ser::AbstractRasterSeries, args...; kw...)
    map(x -> resample(x, args...; kw...), ser)
end
function resample(xs::Union{Tuple,NamedTuple}; to=first(xs), kw...)
    map(x -> resample(x; to, kw...), xs)
end
function resample(A::RasterStackOrArray, resolution::Number; kw...)
    dimres = map(d -> DD.basetypeof(d)(resolution), dims(A, (XDim, YDim)))
    resample(A, dimres; kw...)
end
function resample(A::RasterStackOrArray, pairs::Pair...; kw...)
    resample(A, map((D, x) -> D(x), pairs); kw...)
end
function resample(A::RasterStackOrArray, resdims::DimTuple;
    crs::GeoFormat=crs(A), method=:near, kw...
)
    wkt = convert(String, convert(WellKnownText, crs))
    res = map(val, dims(resdims, (YDim, XDim)))
    flags = Dict(
        :t_srs => wkt,
        :tr => res,
        :r => method,
    )
    return warp(A, flags; kw...)
end
function resample(A::RasterStackOrArray; to, method=:near, kw...)
    all(hasdim(to, (XDim, YDim))) || throw(ArgumentError("`to` mush have both XDim and YDim dimensions to resize with GDAL"))
    if sampling(to, XDim) isa Points
        to = set(to, dims(to, XDim) => Intervals(Start()))
    end
    if sampling(to, YDim) isa Points
        to = set(to, dims(to, YDim) => Intervals(Start()))
    end

    wkt = convert(String, convert(WellKnownText, crs(to)))
    xres, yres = map(abs ∘ step, span(to, (XDim, YDim)))
    (xmin, xmax), (ymin, ymax) = bounds(to, (XDim, YDim))
    flags = Dict(
        :t_srs => wkt,
        :tr => [yres, xres],
        :te => [xmin, ymin, xmax, ymax],
        :r => method,
    )
    return warp(A, flags; kw...)
end
