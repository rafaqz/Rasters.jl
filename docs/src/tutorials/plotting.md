## Plots, simple
[`Plots.jl`](https://github.com/JuliaPlots/Plots.jl) and [`Makie.jl`](https://makie.org) are fully supported by
Rasters.jl, with recipes for plotting `Raster` and `RasterStack` provided. `plot`
will plot a heatmap with axes matching dimension values. If `mappedcrs` is used,
converted values will be shown on axes instead of the underlying `crs` values.
`contourf` will similarly plot a filled contour plot.

Pixel resolution is limited to allow loading very large files quickly. `max_res` 
specifies the maximum pixel resolution to show on the longest axis of the array.
It can be set manually to change the resolution (e.g. for large or high-quality plots):

````@example plots
using Rasters, RasterDataSources, ArchGDAL, Plots
A = Raster(WorldClim{BioClim}, 5)
plot(A; max_res=3000)
````

For Makie, `plot` functions in a similar way.  `plot` will only accept two-dimensional rasters.
You can invoke `contour`, `contourf`, `heatmap`, `surface` or any Makie plotting function which
supports surface-like data on a **2D raster**.

To obtain tiled plots for 3D rasters and RasterStacks, use the function `Rasters.rplot([gridposition], raster; kw_args...)`. 
This is an unexported function, since we're not sure how the API will change going forward.

## Makie, simple

````@example plots
using CairoMakie 
CairoMakie.activate!(px_per_unit = 2)
using Rasters, CairoMakie, RasterDataSources, ArchGDAL
A = Raster(WorldClim{BioClim}, 5)
Makie.plot(A)
````

## Loading data
Our first example simply loads a file from disk and plots it.

This `netcdf` file only has one layer, if it has more we could use RasterStack instead.

````@example plots
using Rasters, NCDatasets, Plots
using Downloads: download

url = "https://archive.unidata.ucar.edu/software/netcdf/examples/tos_O1_2001-2002.nc";
filename = download(url, "tos_O1_2001-2002.nc");
A = Raster(filename)
````

Objects with Dimensions other than X and Y will produce multi-pane plots. Here we plot every third month in
the first year in one plot:

````@example plots
A[Ti=1:3:12] |> plot
````

Now plot the ocean temperatures around the Americas in the first month of 2001.
Notice we are using lat/lon coordinates and date/time instead of regular
indices. The time dimension uses `DateTime360Day`, so we need to load CFTime.jl
to index it with `Near`.

````@example plots
using CFTime
A[Ti(Near(DateTime360Day(2001, 01, 17))), Y(-60.0 .. 90.0), X(45.0 .. 190.0)] |> plot 
````

Now get the mean over the timespan, then save it to disk, and plot it as a
filled contour.

Other plot functions and sliced objects that have only one `X`/`Y`/`Z` dimension
fall back to generic DimensionalData.jl plotting, which will still correctly
label plot axes.

````@example plots
using Statistics
# Take the mean
mean_tos = mean(A; dims=Ti)
````

### Plot a contour plot

````@example plots
using Plots
Plots.contourf(mean_tos; dpi=300, size=(800, 400))
````
### write to disk
Write the mean values to disk

````@example plots
write("mean_tos.nc", mean_tos)
````

Plotting recipes in DimensionalData.jl are the fallback for Rasters.jl when the
object doesn't have 2 `X`/`Y`/`Z` dimensions, or a non-spatial plot command is
used. So (as a random example) we could plot a transect of ocean surface
temperature at 20 degree latitude :

````@example plots
A[Y(Near(20.0)), Ti(1)] |> plot
````
## Polygon masking, mosaic and plot

In this example we will `mask` the Scandinavian countries with border polygons,
then `mosaic` together to make a single plot. 

First, get the country boundary shape files using GADM.jl.

using Rasters, RasterDataSources, ArchGDAL, Shapefile, Plots, Dates, Downloads, NCDatasets

## Download the shapefile
````@example plots
using Downloads
using Shapefile
shapefile_url = "https://github.com/nvkelso/natural-earth-vector/raw/master/10m_cultural/ne_10m_admin_0_countries.shp"
shapefile_name = "boundary_lines.shp"
Downloads.download(shapefile_url, shapefile_name);
nothing # hide
````
### Load using Shapefile.jl

````@example plots
shapes = Shapefile.Handle(shapefile_name)
denmark_border = shapes.shapes[71]
norway_border = shapes.shapes[53]
sweden_border = shapes.shapes[54];
nothing # hide
````

Then load raster data. We load some worldclim layers using `RasterDataSources` via Rasters.jl:

````@example plots
using Rasters, RasterDataSources
using Dates
climate = RasterStack(WorldClim{Climate}, (:tmin, :tmax, :prec, :wind); month=July)
````

`mask` Denmark, Norway and Sweden from the global dataset using their border polygon,
then trim the missing values. We pad `trim` with a 10 pixel margin.

````@example plots
mask_trim(climate, poly) = trim(mask(climate; with=poly); pad=10)

denmark = mask_trim(climate, denmark_border)
norway = mask_trim(climate, norway_border)
sweden = mask_trim(climate, sweden_border)
````

### Plotting with Plots.jl
First define a function to add borders to all subplots.

````@example plots
function borders!(p, poly)
    for i in 1:length(p)
        Plots.plot!(p, poly; subplot=i, fillalpha=0, linewidth=0.6)
    end
    return p
end
````

Now we can plot the individual countries.

````@example plots
dp = plot(denmark)
borders!(dp, denmark_border)
````

and sweden

````@example plots
sp = plot(sweden)
borders!(sp, sweden_border)
````

and norway

````@example plots
np = plot(norway)
borders!(np, norway_border)
````

The Norway shape includes a lot of islands. Lets crop them out using `..` intervals:

````@example plots
norway_region = climate[X(0..40), Y(55..73)]
plot(norway_region)
````
And mask it with the border again:

````@example plots
norway = mask_trim(norway_region, norway_border)
np = plot(norway)
borders!(np, norway_border)
````

Now we can combine the countries into a single raster using mosaic. first will take the first value if/when there is an overlap.

````@example plots
scandinavia = mosaic(first, denmark, norway, sweden)
````

And plot scandinavia, with all borders included:

````@example plots
p = plot(scandinavia)
borders!(p, denmark_border)
borders!(p, norway_border)
borders!(p, sweden_border)
p
````

And save to `netcdf` - a single multi-layered file, and `tif`, which will write a file for each stack layer.

````@example plots
write("scandinavia.nc", scandinavia)
write("scandinavia.tif", scandinavia)
````

`Rasters.jl` provides a range of other methods that are being added to over time. Where applicable these methods
read and write lazily to and from disk-based arrays of common raster file types. These methods also work for
entire RasterStacks and RasterSeries using the same syntax.
