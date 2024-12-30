```@meta
CollapsedDocStrings=true
```
# Reprojection and resampling

### What is resampling?

**[`resample`](@ref)** "re-samples" the 
data by interpolation and can also aggregate or disaggregate, changing the resolution.
It always returns a `Regular` lookup (like a range), and is the most flexible of the 
resampling methods.

This uses GDAL's `gdalwarp` algorithm under the hood.  You can call that via [`warp`](@ref)
if you need more control, but generally `resample` is sufficient. 

::: tip warp, contributions are welcome!

- Show how to use `warp` to reproject a raster

:::

Rasters.jl has a few other methods to change the lookups of a raster.  These are:
- [`reproject`](@ref), which directly reprojects the lookup axes 
  (but is **only usable for specific cases**, where the source and destination 
  coordinate systems are both cylindrical, like the long-lat, Mercator, or Web-Mercator projections.) 

  This is a lossless operation and keeps the data exactly the same - only the axes are changed. 

- [`aggregate`](@ref) and [`disaggregate`](@ref), which change the resolution of 
  the raster by merging ([`aggregate`](@ref)) or splitting ([`disaggregate`](@ref)) cells.

  They can't change cells fractionally, and can't change the projection or coordinate system.

Of all these methods, **`resample`** is the most flexible and powerful, and is what we will focus on here.  
It is, however, also the slowest.  So if another method is applicable to your problem, you should consider it.

### How `resample` works

`resample` uses GDAL's `gdalwarp` algorithm under the hood.  This is a battle-tested algorithm
and is generally pretty robust.  However, it has the following limitations:
- It always assumes cell-based sampling, instead of point-based sampling.  This does mean that 
  point-based rasters are converted to cell-based sampling.
- It can only accept some primitive types for the input data, since that data is passed directly to a C library.
  Things like `RGB` or user-defined types are not usually supported.

`resample` allows you to specify several methods, see some of them in the next section.

### `resolution`, `size` and `methods`

Let's start by loading the necessary packages:

````@example resample
using Rasters, RasterDataSources, ArchGDAL
using DimensionalData
using DimensionalData.Lookups
using NaNStatistics
using CairoMakie
````

````@example resample
ras = Raster(WorldClim{BioClim}, 5)
ras_m = replace_missing(ras, missingval=NaN)
````

resampling to a given `size` or `res ≡ resolution` providing a `method` is done with

:::tabs

== tab size

````@ansi resample
ras_sample = resample(ras_m; size=(2160, 1080), method="average")
````

== tab resolution

````@ansi resample
ras_sample = resample(ras_m; res=1.0, method="average")
````

:::

other available methods to try:  `"mode"`, `"max"`, `"sum"`, `"bilinear"`, `"cubic"`, `"cubicspline"`, `"lanczos"`, `"min"`, `"med"`, `"q1"`, `"q3"` and  `"near"`.

Let's consider a few more examples, with the following options:

````@example resample
methods = ["average", "mode", "max", "sum"]
sizes = [(2160, 1080), (1440, 720), (720, 360), (360, 180)]
resolutions = [0.16666666666666666, 0.25, 0.5, 1.0];
nothing # hide
````

:::tabs

== tab sizes and methods

````@example resample
method_sizes = [resample(ras_m; size=size, method=method) for method in methods for size in sizes]
with_theme(Rasters.theme_rasters()) do
    colorrange = (nanminimum(ras_m), nanmaximum(ras_m))
    hm=nothing
    fig = Figure(; size = (1000, 600))
    axs = [Axis(fig[i,j], title="size=$(size), method=:$(method)", titlefont=:regular)
        for (i, method) in enumerate(methods) for (j, size) in enumerate(sizes)]
    for (i, ax) in enumerate(axs)
        hm = heatmap!(ax, method_sizes[i]; colorrange)
    end
    Colorbar(fig[:,end+1], hm)
    hidedecorations!.(axs; grid=false)
    rowgap!(fig.layout, 5)
    colgap!(fig.layout, 10)
    fig
end
````

== tab resolutions and methods

````@example resample
method_res = [resample(ras_m; res=res, method=method) for method in methods for res in resolutions]
with_theme(Rasters.theme_rasters()) do
    colorrange = (nanminimum(ras_m), nanmaximum(ras_m))
    hm=nothing
    fig = Figure(; size = (1000, 600))
    axs = [Axis(fig[i,j], title="res=$(round(res, digits=4)), method=:$(method)", titlefont=:regular)
        for (i, method) in enumerate(methods) for (j, res) in enumerate(resolutions)]
    for (i, ax) in enumerate(axs)
        hm = heatmap!(ax, method_res[i]; colorrange)
    end
    Colorbar(fig[:,end+1], hm)
    hidedecorations!.(axs; grid=false)
    rowgap!(fig.layout, 5)
    colgap!(fig.layout, 10)
    fig
end
````

:::


## `reproject` with `resample` using a `ProjString`

Geospatial datasets will come in different [projections](https://proj.org/en/9.4/operations/projections/index.html) or coordinate reference systems (CRS) for many reasons. Here, we will focus on `MODIS SINUSOIDAL` and `EPSG`, and transformations between them.

Let's load our test raster

````@example resample
ras = Raster(WorldClim{BioClim}, 5)
ras_m = replace_missing(ras, missingval=NaN);
nothing # hide
````

### Sinusoidal Projection (MODIS)

````@example resample
SINUSOIDAL_CRS = ProjString("+proj=sinu +lon_0=0 +type=crs")
````

::: details Raw MODIS ProjString

````julia
SINUSOIDAL_CRS = ProjString("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
````

:::


and the `resample` is performed with

````@example resample
ras_sin = resample(ras_m; size=(2160, 1080), crs=SINUSOIDAL_CRS, method="average") # ? do again, failing locally on first try
ras_sin = resample(ras_m; size=(2160, 1080), crs=SINUSOIDAL_CRS, method="average") # hide
nothing # hide
````

::: tip

`GDAL` always changes the locus to cell sampling, you can reset this by using `shiftlocus`.

:::

let's compare the total counts!

````@example resample
nansum(ras_m), nansum(ras_sin)
````

and, how does this looks like?

````@example resample
fig, ax, plt = heatmap(ras_sin)
Colorbar(fig[1,2], plt)
fig
````

now, let's go back to `latitude` and `longitude` and reduce the resolution

````@example resample
ras_epsg = resample(ras_sin; size=(1440,720), crs=EPSG(4326), method="average")
````

and let's apply `shiftlocus` such that the lookups share the exact same grid, which might be needed when building bigger datasets:

````@example resample
locus_resampled = DimensionalData.shiftlocus(Center(), ras_epsg)
````

and compare the total counts!

````@example resample
nansum(ras_m), nansum(locus_resampled)
````

::: info Things to keep in mind

  - You can in fact resample to another raster `resample(ras; to=ref_ras)`, if you want perfect alignment. Contributions are welcome for this use case!
  - This doesn't work for irregularly sampled rasters.

:::


````@example resample
fig, ax, plt = heatmap(ras_epsg)
Colorbar(fig[1,2], plt)
fig
````

### A `Raster` from scratch

````@example resample
x_range = LinRange(-180, 179.75, 1440)
y_range = LinRange(89.75, -90, 720)
ras_data = ras_epsg.data
nothing # hide
````

create the raster

````@ansi modis
ras_scratch = Raster(ras_data, (X(x_range; sampling=Intervals(Start())),
    Y(y_range; sampling=Intervals(Start()))), crs=EPSG(4326), missingval=NaN)

````

::: warning

Note that you need to specify `sampling=Intervals(Start())` for `X` and `Y`.

This requires that you run `using Rasters.Lookups`, where the `Intervals` and `Start` types are defined.

:::

and take a look

````@example resample
fig, ax, plt = heatmap(ras_scratch)
Colorbar(fig[1,2], plt)
fig
````

and the corresponding resampled projection

````@ansi modis
ras_sin_s = resample(ras_scratch; size=(1440,720), crs=SINUSOIDAL_CRS, method="average")
````

````@example resample
fig, ax, plt = heatmap(ras_sin_s)
Colorbar(fig[1,2], plt)
fig
````

and go back from `sin` to `epsg`:

````@example resample
ras_epsg = resample(ras_sin_s; size=(1440,720), crs=EPSG(4326), method="average")
locus_resampled = DimensionalData.shiftlocus(Center(), ras_epsg)

fig, ax, plt = heatmap(locus_resampled)
Colorbar(fig[1,2], plt)
fig
````


and compare the total counts again!

````@example resample
nansum(ras_m), nansum(locus_resampled)
````

::: danger

Note that all counts are a little bit off. Could we mitigate this some more?

:::



