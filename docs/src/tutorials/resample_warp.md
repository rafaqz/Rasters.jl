```@meta
CollapsedDocStrings=true
```
# Reprojection and resampling

- What is resampling?
    - When to resample vs reproject
    - Things to keep in mind
        - GDAL always changes the locus to cell sampling, you can reset this by using `shiftlocus`
        - You can in fact resample to another raster, if you want perfect alignment.
            - This doesn't work for irregularly sampled rasters.
        - Resampling is a lossy operation and takes time.  Try to avoid repeatedly resampling, and if you must do so, crop or trim the raster as much as you can first.
- Show the different resampling methods, in a grid
- Show some different projections and ways of constructing them
- Show how to use `size` and `res` to change the resolution of a raster
- Show how to use `warp` to reproject a raster

## What is resampling?

**[`resample`](@ref)** "re-samples" the 
data by interpolation and can also aggregate or disaggregate, changing the resolution.
It always returns a `Regular` lookup (like a range), and is the most flexible of the 
resampling methods.

This uses GDAL's `gdalwarp` algorithm under the hood.  You can call that via [`warp`](@ref)
if you need more control, but generally `resample` is sufficient. 

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

## How `resample` works

`resample` uses GDAL's `gdalwarp` algorithm under the hood.  This is a battle-tested algorithm
and is generally pretty robust.  However, it has the following limitations:
- It always assumes cell-based sampling, instead of point-based sampling.  This does mean that 
  point-based rasters are converted to cell-based sampling.
- It can only accept some primitive types for the input data, since that data is passed directly to a C library.
  Things like `RGB` or user-defined types are not usually supported.

`resample` allows you to specify several methods, which you can see if you expand the docstring below.

## `resample` and projections with `ProjString`

Geospatial datasets will come in different [projections](https://proj.org/en/9.4/operations/projections/index.html) or coordinate reference systems (CRS) for many reasons. Here, we will focus on `MODIS SINUSOIDAL` and `EPSG`, and transformations between them.

Let's start by loading the neccesary packages:

````@example modis
using Rasters, RasterDataSources, ArchGDAL
using DimensionalData
using DimensionalData.Lookups
using NaNStatistics
using CairoMakie
````

and let's load a test raster

````@example modis
ras = Raster(WorldClim{BioClim}, 5)
ras_m = replace_missing(ras, missingval=NaN)
````

and let's also take a look

````@example modis
fig, ax, plt = heatmap(ras_m)
Colorbar(fig[1,2], plt)
fig

````

### MODIS SINUSOIDAL PROJECTION

````@example modis
# ? is this the right ProjString ?
# SINUSOIDAL_CRS = ProjString("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
SINUSOIDAL_CRS = ProjString("+proj=sinu +lon_0=0 +type=crs")
````

and hence the `resample` is performed with

````@example modis
ras_sin = resample(ras_m; size=(2160, 1080), crs=SINUSOIDAL_CRS, method="average") # ? do again, failing locally on first try
ras_sin = resample(ras_m; size=(2160, 1080), crs=SINUSOIDAL_CRS, method="average") # hide
nothing # hide
````

let's compare the total counts!

````@example modis
nansum(ras_m), nansum(ras_sin)
````

and, how does this looks like?

````@example modis
fig, ax, plt = heatmap(ras_sin)
Colorbar(fig[1,2], plt)
fig
````

now, let's go back to `latitude` and `longitude` and reduce the resolution

````@example modis
ras_epsg = resample(ras_sin; size=(1440,720), crs=EPSG(4326), method="average")
````

and let's apply `shiftlocus` such that the lookups share the exact same grid, which might be needed when building bigger datasets:

````@example modis
locus_resampled = DimensionalData.shiftlocus(Center(), ras_epsg)
````

and compare the total counts!

````@example modis
nansum(ras_m), nansum(locus_resampled)
````

````@example modis
fig, ax, plt = heatmap(ras_epsg)
Colorbar(fig[1,2], plt)
fig
````

### Construct a Raster from scratch natively in the sinusoidal projection

````@example modis
x_range = LinRange(-180, 179.75, 1440)
y_range = LinRange(89.75, -90, 720)
ras_data = ras_epsg.data
nothing # hide
````

create the raster

````@ansi modis
ras_scratch = Raster(ras_data, (X(x_range; sampling=Intervals(Start())),
    Y(y_range; sampling=Intervals(Start()))), crs=EPSG(4326))

````

::: warning

Note that you need to specify `sampling=Intervals(Start())` for `X` and `Y`.

This requires that you run `using Rasters.Lookups`, where the `Intervals` and `Start` types are defined.

:::

and take a look

````@example modis
fig, ax, plt = heatmap(ras_scratch)
Colorbar(fig[1,2], plt)
fig
````

and the corresponding resampled projection

````@ansi modis
ras_sin_s = resample(ras_scratch; size=(1440,720), crs=SINUSOIDAL_CRS, method="average")
````

````@example modis
fig, ax, plt = heatmap(ras_sin_s)
Colorbar(fig[1,2], plt)
fig
````

and go back from `sin` to `epsg`:

````@example modis
ras_epsg = resample(ras_sin_s; size=(1440,720), crs=EPSG(4326), method="average")
locus_resampled = DimensionalData.shiftlocus(Center(), ras_epsg)

fig, ax, plt = heatmap(locus_resampled)
Colorbar(fig[1,2], plt)
fig
````


and compare the total counts again!

````@example modis
nansum(ras_m), nansum(locus_resampled)
````

::: danger

Note that all counts are a little bit off. Could we mitigate this some more?

:::



