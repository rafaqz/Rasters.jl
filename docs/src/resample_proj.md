## resample and Projections with ProjString

Geospatial datasets will come in different [projections](https://proj.org/en/9.4/operations/projections/index.html) for whaever reasons. Here, we will focus in some of the most used ones, e.g., `MODIS SINUSOIDAL` and `EPSG`, as well as transformations between them.

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
fig = plot(ras_m; colorrange=(0,100))
````

### MODIS SINUSOIDAL PROJECTION

````@example
# ? is this the right ProjString ?, do we need to shift lat, lon?
SINUSOIDAL_CRS = ProjString("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")
````

and hence the `resample` is performed with

````@example modis
ras_sin = resample(ras_m; size=(1440,720), crs=SINUSOIDAL_CRS, method="sum")
````

let's compare the total counts!

````@example modis
nansum(ras_m), nansum(ras_sin)
````

and, how does this looks like?

````@example modis
fig = plot(ras_sin; colorrange=(0,100))
````

now, let's go back to `latitude` and `longitude`

````@example modis
# ? also here, do we need to shift X, Y?
ras_epsg = resample(ras_sin; size=(1440,720), crs=EPSG(4326), method="sum")
````

and let's apply `shiftlocus` such that we can harmonize coordinates, which might be needed when building bigger datasets:

````@example modis
locus_resampled = DimensionalData.shiftlocus(Center(), ras_epsg)
````

and compare the total counts!

````@example modis
nansum(ras_m), nansum(locus_resampled)
````

````@example modis
fig = plot(ras_epsg; colorrange=(0,100))
````

### Assemble a Raster from scratch nativily in the sinusoidal projection

````@example modis
x_range = range(-2.0015109355797417e7, 1.998725401355172e7, 1440)
y_range = range(9.979756529777847e6, -1.0007111969122082e7, 720)
ra_data = ras_sin.data;
nothing # hide
````

create the raster

````@example modis
ras_scratch = Raster(ra_data, (X(x_range; sampling=Intervals(Start())), Y(y_range; sampling=Intervals(Start()))),
    crs=SINUSOIDAL_CRS)
````

::: warning
At the moment, you need to specify `sampling=Intervals(Start())` for `X` and `Y`.
:::

and take a look

````@example modis
fig = plot(ras_scratch; colorrange=(0,100))
````

and the corresponding resampled projection

````@example modis
ras_latlon = resample(ras_scratch; size=(1440,720), crs=EPSG(4326), method="sum")
locus_resampled  = DimensionalData.shiftlocus(Center(), ras_latlon)
````

````@example modis
fig = plot(ras_latlon; colorrange=(0,100))
````

and compare the total counts again!

````@example modis
nansum(ras_m), nansum(locus_resampled)
````

::: danger
Note that all counts are a little bit off. Could we mitigate this a some more?
:::



