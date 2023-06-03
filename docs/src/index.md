# Rasters.jl

```@docs
Rasters
```

# Common Applications

## Subsetting an object

Regular `getindex` (e.g. `A[1:100, :]`) and `view` work on all objects just as
with an `Array`. `view` is always lazy, and reads from disk are deferred until
`getindex` is used. DimensionalData.jl `Dimension`s and `Selector`s are the other
way to subset an object, making use of the objects index to find values at 
e.g. certain X/Y coordinates. The available selectors are listed here:

|                        |                                                                    |
| :--------------------- | :----------------------------------------------------------------- |
| `At(x)`                | get the index exactly matching the passed in value(s).             |
| `Near(x)`              | get the closest index to the passed in value(s).                   |
| `Where(f::Function)`   | filter the array axis by a function of the dimension index values. |
| `a..b`/`Between(a, b)` | get all indices between two values, excluding the high value.      |
| `Contains(x)`          | get indices where the value x falls within an interval.            |


Use the `..` selector to take a `view` of madagascar:

```@example
using Rasters, RasterDataSources, ArchGDAL, Plots
A = Raster(WorldClim{BioClim}, 5)
madagascar = view(A, X(43.25 .. 50.48), Y(-25.61 .. -12.04)) # Note the space between .. -12
plot(madagascar)
```

## Methods that change the reslolution or extent of an object

Click through to the function documentation for more in-depth descriptions and
examples.

|                           |                                                                              |
| :------------------------ | :--------------------------------------------------------------------------- |
| [`aggregate`](@ref)       | aggregate data by the same or different amounts for each axis.               |
| [`disaggregate`](@ref)    | similarly disaggregate data.                                                 |
| [`mosaic`](@ref)          | join rasters covering different extents into a single array or file.         |
| [`crop`](@ref)            | shrink objects to specific dimension sizes or the extent of another object.  |
| [`extend`](@ref)          | extend objects to specific dimension sizes or the extent of another object.  |
| [`trim`](@ref)            | trims areas of missing values for arrays and across stack layers.            |
| [`resample`](@ref)        | resample data to a different size and projection, or snap to another object. |
| [`warp`](@ref)            | use `gdalwarp` on any object, e.g. a multidimensional NetCDF stack.          |


## Methods that change an objects values: 

Note that most regular Julia methods, such as `replace`, work as for a standard
`Array`. These additional methods are commonly required in GIS applications.

|                           |                                                                              |
| :------------------------ | :--------------------------------------------------------------------------- |
| [`classify`](@ref)        | classify values into categories.                                             |
| [`mask`](@ref)            | mask an object by a polygon or `Raster` along `X/Y`, or other dimensions.    |
| [`replace_missing`](@ref) | replace all missing values in an object and update `missingval`.             |


## Point, polygon and table operation

|                           |                                                                              |
| :------------------------ | :--------------------------------------------------------------------------- |
| [`rasterize`](@ref)       | rasterize points and geometries.                                             |
| [`coverage`](@ref)        | get the fraction of each pixel covered by geometries.                        |
| [`extract`](@ref)         | extract values from points or geometries.                                    |
| [`zonal`](@ref)           | calculate zonal statistics for an object masked by geometries.               |


## Methods to load, write and modify data sources:

|                           |                                                                         |
| :------------------------ | :---------------------------------------------------------------------- |
| [`modify`](@ref)          | replace the data in objects. Useful to e.g. move objects to/from a GPU. |
| [`read`](@ref)            | read data to memory if it is on disk.                                   |
| [`read!`](@ref)           | read data to predefined memory.                                         |
| [`open`](@ref)            | open the underlying data for manually reading or writing.               |
| [`write`](@ref)           | write objects to file.                                                  |


## Altering and summarising arrays and stacks with regular julia methods

Most base methods work as for regular julia `Array`s, such as `reverse` and
rotations like `rotl90`. Base, statistics and linear algebra methods like `mean`
that take a `dims` argument can also use the dimension name. To take the mean
over the time dimension:

```julia
mean(A, dims=Ti)
```

`broadcast` works lazily from disk when `lazy=true`, and is only applied when data
is directly indexed. Adding a dot to any function will use broadcast over a `Raster`
just like an `Array`. 

## Broadcasting

For a disk-based array `A`, this will only be applied when indexing occurs or
when we [`read`](@ref) the array.

```julia
A .*= 2
```

To broadcast directly to disk, we need to open the file in write mode:

```julia
open(Raster(filename); write=true) do O
    O .*= 2
end
```

To broadcast over a `RasterStack` use `map`, which applies a function to
the raster layers of the stack.

```julia
newstack = map(stack) do raster
    raster .* 2
end
```

## Modifying object properties

`rebuild` can be used to modify the fields of an object, generating a new object
(but possibly holding the same arrays or files).

If you know that a file had an incorrectly specified missing value, you can do:

```julia
rebuild(A; missingval=-9999)
```

(`replace_missing` will actualy replace the current values)

Or if you need to change the name of the layer:

```julia
rebuild(A; name=:temperature)
```

`set` can be used to modify the nested properties of an objects dimensions, that
are more difficult to change with `rebuild`. `set` works on the principal that
dimension properties can only be in one specific field, so we generally don't
have to specify which one it is. `set` will also try to update anything affected
by a change you make.

This will set the `X` axis to specify points, instead of intervals:

```julia
set(A, X => Points)
```

We can also reassign dimensions, here `X` becomes `Z`:

```julia
set(A, X => Z)
```

`setcrs(A, crs)` and `setmappedcrs(A, crs)` will set the crs value/s of an
object to any `GeoFormat` from GeoFormatTypes.jl.


# Examples and Plotting

[Plots.jl](https://github.com/JuliaPlots/Plots.jl) and [Makie.jl](https://github.com/MakieOrg/Plots.jl) are fully supported by
Rasters.jl, with recipes for plotting `Raster` and `RasterStack` provided. `plot`
will plot a heatmap with axes matching dimension values. If `mappedcrs` is used,
converted values will be shown on axes instead of the underlying `crs` values.
`contourf` will similarly plot a filled contour plot.

Pixel resolution is limited to allow loading very large files quickly. `max_res` 
specifies the maximum pixel resolution to show on the longest axis of the array.
It can be set manually to change the resolution (e.g. for large or high-quality plots):

```julia
using Rasters, RasterDataSources, ArchGDAL, Plots
A = Raster(WorldClim{BioClim}, 5)
plot(A; max_res=3000)
```

For Makie, `plot` functions in a similar way.  `plot` will only accept two-dimensional rasters.  You can invoke `contour`, `contourf`, `heatmap`, `surface` or any Makie plotting function which supports surface-like data on a **2D raster**.


To obtain tiled plots for 3D rasters and RasterStacks, use the function `Rasters.rplot([gridposition], raster; kw_args...)`.  This is an unexported function, since we're not sure how the API will change going forward.

```@example makie
using CairoMakie # hide
CairoMakie.activate!(px_per_unit = 2) # hide
using Rasters, CairoMakie, RasterDataSources, ArchGDAL
A = Raster(WorldClim{BioClim}, 5)
Makie.plot(A)
```

## Loading and plotting data

Our first example simply loads a file from disk and plots it.

This netcdf file only has one layer, if it has more we could use `RasterStack`
instead. 

```@example nc
using Rasters, NCDatasets, Plots
url = "https://www.unidata.ucar.edu/software/netcdf/examples/tos_O1_2001-2002.nc";
filename = download(url, "tos_O1_2001-2002.nc");
A = Raster(filename)
```

Objects with Dimensions other than `X` and `Y` will produce multi-pane plots.
Here we plot every third month in the first year in one plot:

```@example nc
A[Ti=1:3:12] |> plot
```

Now plot the ocean temperatures around the Americas in the first month of 2001.
Notice we are using lat/lon coordinates and date/time instead of regular
indices. The time dimension uses `DateTime360Day`, so we need to load CFTime.jl
to index it with `Near`.

```@example nc
using CFTime
A[Ti(Near(DateTime360Day(2001, 01, 17))), Y(-60.0 .. 90.0), X(45.0 .. 190.0)] |> plot 
```

Now get the mean over the timespan, then save it to disk, and plot it as a
filled contour:

Other plot functions and sliced objects that have only one `X`/`Y`/`Z` dimension
fall back to generic DimensionalData.jl plotting, which will still correctly
label plot axes.

```@example nc
using Statistics
# Take the mean
mean_tos = mean(A; dims=Ti)
```

## Plot a contour plot

```@example nc
contourf(mean_tos; dpi=300, size=(800, 400))
```

Write the mean values to disk

```@example nc
write("mean_tos.nc", mean_tos)
```

Plotting recipes in DimensionalData.jl are the fallback for Rasters.jl when the
object doesn't have 2 `X`/`Y`/`Z` dimensions, or a non-spatial plot command is
used. So (as a random example) we could plot a transect of ocean surface
temperature at 20 degree latitude :

```@example nc
A[Y(Near(20.0)), Ti(1)] |> plot
```

## A basic species distribution modelling workflow

Load occurrences for the Mountain Pygmy Possum using GBIF.jl

```@example sdm
using Rasters, RasterDataSources, ArchGDAL, GBIF2, Plots 
records = GBIF2.occurrence_search("Burramys parvus"; limit=300)
```

Extract the longitude/latitude value to a `Vector` of points
(a `Tuple` counts as a `(x, y)` point in GeoInterface.jl):

```@example sdm
coords = [(r.decimalLongitude, r.decimalLatitude) for r in records if !ismissing(r.decimalLatitude)]
```

Get BioClim layers and subset to south-east Australia

```@example sdm
A = RasterStack(WorldClim{BioClim}, (1, 3, 7, 12))
se_aus = A[X(138 .. 155), Y(-40 .. -25)]
```

Plot BioClim predictors and scatter occurrence points on all subplots

```@example sdm
p = plot(se_aus);
kw = (legend=:none, opacity=0.5, markershape=:cross, markercolor=:black)
foreach(i -> scatter!(p, coords; subplot=i, kw...), 1:4)
display(p)
```

Then extract predictor variables and write to CSV.

```@example sdm
using CSV
predictors = collect(extract(se_aus, coords))
CSV.write("burramys_parvus_predictors.csv", predictors)
```

Or convert them to a `DataFrame`.

```@example sdm
using DataFrames
df = DataFrame(predictors)
df[1:5, :]
```

## Polygon masking, mosaic and plot

In this example we will `mask` the Scandinavian countries with border polygons,
then `mosaic` together to make a single plot. 

First, get the country boundary shape files using GADM.jl.

```@example mask
using Rasters, RasterDataSources, ArchGDAL, Shapefile, Plots, Dates, Downloads, NCDatasets

# Download the shapefile
shapefile_url = "https://github.com/nvkelso/natural-earth-vector/raw/master/10m_cultural/ne_10m_admin_0_countries.shp"
shapefile_name = "boundary_lines.shp"
Downloads.download(shapefile_url, shapefile_name)

# Load using Shapefile.jl
shapes = Shapefile.Handle(shapefile_name)
denmark_border = shapes.shapes[71]
norway_border = shapes.shapes[53]
sweden_border = shapes.shapes[54]
```

Then load raster data. We load some worldclim layers using `RasterDataSources` via Rasters.jl:

```@example mask
climate = RasterStack(WorldClim{Climate}, (:tmin, :tmax, :prec, :wind); month=July)
```

`mask` Denmark, Norway and Sweden from the global dataset using their border polygon,
then trim the missing values. We pad `trim` with a 10 pixel margin.

```@example mask
mask_trim(climate, poly) = trim(mask(climate; with=poly); pad=10)

denmark = mask_trim(climate, denmark_border)
norway = mask_trim(climate, norway_border)
sweden = mask_trim(climate, sweden_border)
```

## Plotting

First define a function to add borders to all subplots.

```@example mask
function borders!(p, poly) 
    for i in 1:length(p)
        plot!(p, poly; subplot=i, fillalpha=0, linewidth=0.6)
    end
    return p
end
```

Now we can plot the individual countries.

```@example mask
dp = plot(denmark)
borders!(dp, denmark_border)
```

```@example mask
sp = plot(sweden)
borders!(sp, sweden_border)
```

```@example mask
np = plot(norway)
borders!(np, norway_border)
```

The Norway shape includes a lot of islands. Lets crop them out using `..` intervals:

```@example mask
norway_region = climate[X(0..40), Y(55..73)]
plot(norway_region)
```

And mask it with the border again:
```@example mask
norway = mask_trim(norway_region, norway_border)
np = plot(norway)
borders!(np, norway_border)
```

Now we can combine the countries into a single raster using `mosaic`.
`first` will take the first value if/when there is an overlap.

```@example mask
scandinavia = mosaic(first, denmark, norway, sweden)
```

And plot scandinavia, with all borders included:

```@example mask
p = plot(scandinavia)
borders!(p, denmark_border)
borders!(p, norway_border)
borders!(p, sweden_border)
```

And save to netcdf - a single multi-layered file, and tif, which will write a
file for each stack layer.

```@example mask
write("scandinavia.nc", scandinavia)
write("scandinavia.tif", scandinavia)
```

Rasters.jl provides a range of other methods that are being added to over time.
Where applicable these methods read and write lazily to and from disk-based
arrays of common raster file types. These methods also work for entire
`RasterStacks` and `RasterSeries` using the same syntax.

## Plotting in Makie

### 2-D rasters in Makie

Plotting in Makie works somewhat differently than Plots, since the recipe system is different.
You can pass a 2-D raster to any surface-like function (`heatmap`, `contour`, `contourf`, or even `surface` for a 3D plot) with ease.

```@example makie
using CairoMakie, Makie
CairoMakie.activate!(px_per_unit = 2) # hide
using Rasters, RasterDataSources, ArchGDAL
A = Raster(WorldClim{BioClim}, 5) # this is a 3D raster, so is not accepted.
B = A[:, :, 1] # this converts to a 2D raster which Makie accepts!
figure = Figure()
plot(figure[1, 1], B)
contour(figure[1, 2], B)
ax = Axis(figure[2, 1]; aspect = DataAspect())
contourf!(ax, B)
surface(figure[2, 2], B) # even a 3D plot works!
figure
```

### 3-D rasters and RasterStacks in Makie

!!! warning
    This interface is experimental, and unexported for that reason.  It may break at any time!

Just as in Plots, 3D rasters are treated as a series of 2D rasters, which are tiled and plotted.  

You can use `Rasters.rplot` to visualize 3D rasters or RasterStacks in this way.  An example is below:

```@example makie
stack = RasterStack(WorldClim{Climate}; month = 1)
Rasters.rplot(stack; Axis = (aspect = DataAspect(),))
```

You can pass any theming keywords in, which are interpreted by Makie appropriately.

The plots seem a little squished here.  We provide a Makie theme which makes text a little smaller and has some other space-efficient attributes:

```@example makie
CairoMakie.set_theme!(Rasters.theme_rasters())
Rasters.rplot(stack)
```

### Plotting with `Observable`s

`Rasters.rplot` should support Observable input out of the box, but the dimensions of that input
must remain the same - i.e., the element names of a RasterStack must remain the same.

```@example makie
stack_obs = Observable(stack)
fig = Rasters.rplot(stack_obs) # `stack` is the WorldClim climate data for January
record(fig, "rplot.mp4", 1:12; framerate = 3) do i
    stack_obs[] = RasterStack(WorldClim{Climate}; month = i)
end 
```
![](rplot.mp4)

```@eval
using Makie
Makie.set_theme!(Makie.minimal_default)
```

```@docs
Rasters.rplot
```


# Objects

## Raster

Spatial raster data is essentially just an `Array`. But `Raster` wrappers
allow treating them as an array that maintains its spatial index, crs and other
metadata through all transformations. This means they can always be plotted and
written to disk after applying most base Julia methods, and most `broadcast`s.

```@docs
AbstractRaster
Raster
Raster(T::Type{<:RasterDataSources.RasterDataSource}, layer)
```

## RasterStack

Spatial data often comes as a bundle of multiple named arrays, as in netcdf.
`RasterStack` can represent this, or multiple files organised in a similar way.

```@docs
AbstractRasterStack
RasterStack
RasterStack(T::Type{<:RasterDataSources.RasterDataSource})
```

## RasterSeries

A series is a meta-array that holds other files/data that is distributed over
some dimension, often time. These files/data can be `Raster`s or `RasterStack`s.

```@docs
AbstractRasterSeries
RasterSeries
RasterSeries(T::Type{<:RasterDataSources.RasterDataSource})
```

## Dimensions

Rasters uses `X`, `Y`, and `Z` dimensions from DimensionalData.jl to represent
spatial directions like longitude, latitude and the vertical dimension, and
subset data with them. `Ti` is used for time, and `Band` represent bands. Other
dimensions can have arbitrary names, but will be treated generically. See
[DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl/) for more
details on how they work. 

```@docs
Band
```

## Lookup Arrays

These specify properties of the index associated with e.g. the X and Y
dimension. Rasters.jl defines additional lookup arrays: [`Projected`](@ref) to handle
dimensions with projections, and [`Mapped`](@ref) where the projection is mapped to
another projection like `EPSG(4326)`. `Mapped` is largely designed to handle
NetCDF dimensions, especially with `Explicit` spans.

```@docs
Rasters.AbstractProjected
Projected
Mapped
```

# Data sources

Rasters.jl uses a number of backends to load raster data. `Raster`, `RasterStack`
and `RasterSeries` will detect which backend to use for you, automatically.

## GRD

R GRD files can be loaded natively, using Julias `MMap` - which means they are
very fast, but are not compressed. They are always 3 dimensional, and have `Y`,
`X` and [`Band`](@ref) dimensions.

## NetCDF

NetCDF `.nc` files are loaded using
[NCDatasets.jl](https://github.com/Alexander-Barth/NCDatasets.jl). Layers from
files can be loaded as `Raster("filename.nc"; key=:layername)`. Without `key`
the first layer is used. `RasterStack("filename.nc")` will use all netcdf variables
in the file that are not dimensions as layers. 

NetCDF layers can have arbitrary dimensions. Known, common dimension names are
converted to `X`, `Y` `Z`, and `Ti`, otherwise `Dim{:layername}` is used. Layers
in the same file may also have different dimensions.

NetCDF files still have issues loading directly from disk for some operations.
Using `read(ncstack)` may help.

## GDAL

All files GDAL can access, such as `.tiff` and `.asc` files, can be loaded,
using [ArchGDAL.jl](https://github.com/yeesian/ArchGDAL.jl/issues). These are
generally best loaded as `Raster("filename.tif")`, but can be loaded as
`RasterStack("filename.tif"; layersfrom=Band)`, taking layers from the `Band`
dimension, which is also the default.

## SMAP

The [Soil Moisture Active-Passive](https://smap.jpl.nasa.gov/) dataset provides
global layers of soil moisture, temperature and other related data, in a custom
HDF5 format. Layers are always 2 dimensional, with `Y` and `X` dimensions.

These can be loaded as multi-layered `RasterStack("filename.h5")`. Individual
layers can be loaded as `Raster("filename.h5"; key=:layerkey)`, without `key`
the first layer is used.

```@docs
smapseries
```

## Writing file formats to disk

Files can be written to disk in all formats other than SMAP HDF5 using
`write("filename.ext", A)`. See the docs for [`write`](@ref). They can (with
some caveats) be written to different formats than they were loaded in as,
providing file-type conversion for spatial data.

Some metadata may be lost in formats that store little metadata, or where
metadata conversion has not been completely implemented.

# RasterDataSources.jl integration

[RasterDataSources.jl](https://github.com/EcoJulia/RasterDataSources.jl)
standardises the download of common raster data sources, with a focus on
datasets used in ecology and the environmental sciences. RasterDataSources.jl is
tightly integrated into Rasters.jl, so that datsets and keywords can be used
directly to download and load data as a `Raster`, `RasterStack`, or `RasterSeries`.

```@example
using Rasters, RasterDataSources, ArchGDAL, Plots, Dates
A = Raster(WorldClim{Climate}, :tavg; month=June)
plot(A)
```

See the docs for [`Raster`](@ref), [`RasterStack`](@ref) and [`RasterSeries`](@ref),
and the docs for `RasterDataSources.getraster` for syntax to specify various
data sources.

# Exported functions

Rasters.jl is a direct extension of DimensionalData.jl. See [DimensionalData.jl
docs](https://rafaqz.github.io/DimensionalData.jl/stable/) for the majority of
types and functions that can be used in Rasters.jl.

Functions more specific to geospatial data are included in Rasters.jl, and
listed below.

```@docs
aggregate
aggregate!
boolmask
classify 
classify!
coverage
coverage!
convertlookup
crop
crs
disaggregate
disaggregate!
extend
extract
mappedcrs
mappedbounds
mappedindex
mask
mask!
missingval
missingmask
mosaic
mosaic!
points
rasterize
rasterize!
resample
replace_missing
reproject
setcrs
setmappedcrs
skipmissing
trim
warp
zonal
```

## Slice and combine to and from `RasterSeries`

```@docs
slice
combine
```

## File operations

These `Base` and `DimensionalData` methods have specific Rasters.jl versions:

```@docs
modify
open
read
read!
write
```

## Internals

```@docs
Rasters.FileArray
Rasters.FileStack
Rasters.OpenStack
Rasters.RasterDiskArray
```
