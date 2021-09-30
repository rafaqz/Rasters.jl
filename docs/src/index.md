# GeoData.jl

```@docs
GeoData
```

## Common Applications

### Subsetting an object

Regular `getindex` (e.g. `A[1:100, :]`) and `view` work on all objects just as
with an `Array`. `view` is always lazy, and reads from disk are deferred until
`getindex` is used. DimensionalData.jl `Dimension`s and `Selector`s are the other
way to subset an object, making use of the objects index to find values at 
e.g. certain X/Y coordinates. The available selectors are listed here:

|                      |                                                                    |
| :------------------- | :----------------------------------------------------------------- |
| `At(x)`              | get the index exactly matching the passed in value(s)              |
| `Near(x)`            | get the closest index to the passed in value(s)                    |
| `Where(f::Function)` | filter the array axis by a function of the dimension index values. |
| `Between(a, b)`      | get all indices between two values, excluding the high value.      |
| `Contains(x)`        | get indices where the value x falls within an interval             |


Use the `Between` selector to take a `view` of madagascar:

```@example
using GeoData, Plots
A = GeoArray(WorldClim{BioClim}, 5)
madagascar = view(A, X(Between(43.25, 50.48)), Y(Between(-25.61, -12.04))) 
plot(madagascar)
savefig("madagascar_bio5.png")
```

![Bioclim 5 for madagascar](madagascar_bio5.png)

### Changing an objects resolution or extent

Methods that change the reslolution or extent of an object are listed here.
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
| [`warp`](@ref)            | use `gdalwarp` directly on any object, e.g. a multidimensional NetCDF stack. |



### Methods that change an objects values: 

Note that most regular Julia methods, such as `replace`, work as for a standard
`Array`. These additional methods are commonly required in GIS applications.

|                           |                                                                              |
| :------------------------ | :--------------------------------------------------------------------------- |
| [`classify`](@ref)        | classify values into categories.                                             |
| [`mask`](@ref)            | mask and object by a polygon or `GeoArray` along `X/Y`, or other dimensions. |
| [`replace_missing`](@ref) | replace all missing values in an object and update `missingval`.             |


### Methods to load, write and modify data sources:

|                           |                                                                         |
| :------------------------ | :---------------------------------------------------------------------- |
| [`modify`](@ref)          | replace the data in objects. Useful to e.g. move objects to/from a GPU. |
| [`read`](@ref)            | read data to memory if it is on disk.                                   |
| [`read!`](@ref)           | read data to predefined memory.                                         |
| [`open`](@ref)            | open the underlying data for manually reading or writing.               |
| [`write`](@ref)           | write objects to file.                                                  |


### Altering and summarising arrays and stacks with regular julia methods

Most base methods work as for regular julia `Array`s, such as `reverse` and
rotations like `rotl90`. Base, statistics and linear algebra methods like `mean`
that take a `dims` argument can also use the dimension name. To take the mean
over the time dimension:

```julia
mean(A dims=Ti)
```

`broadcast` works lazily from disk, and is only applied when data is directly
indexed. Adding a dot to any function will use broadcast over a `GeoArray`. 

### Broadcasting

For a disk-based array `A`, this will only be applied when indexing occurs or
when we [`read`](@ref) the array.

```julia
A .*= 2
```

To broadcast directly to disk, we need to open the file in write mode:

```julia
open(GeoArray(filename); write=true) do O
    O .*= 2
end
```

To broadcast over a `GeoStack` use `map`, which applies a function to the layers
of the stack - here `A`.

```julia
newstack = map(stack) do A
    A .* 2
end
```


### Modifying object properties

`rebuild` can be used to modify the fields of an object, generating a new object
(but possibly holding the same arrays or files).

If you know that a file had an incorrectly specified missing value, you can do:

```julia
rebuild(A; missingval=-9999)
```

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

`setcrs(A, crs)` and `setmappedcrs(A, crs)` will set the crs values to any
`GeoFormat` from GeoFormatTypes.jl. These can't be set directly with `set`, as
they can be the same objects.


## Examples and Plotting

[Plots.jl](https://github.com/JuliaPlots/Plots.jl) is fully supported by
GeoData.jl, with recipes for plotting `GeoArray` and `GeoStack` provided. `plot`
will plot a heatmap with axes matching dimension values. If `mappedcrs` is used,
converted values will be shown on axes instead of the underlying `crs` values.
`contourf` will similarly plot a filled contour plot.

Pixel resolution is limited to allow loading very large files quickly. `max_res` 
specifies the maximum pixel resolution to show on the longest axis of the array.
It can be set manually to change the resolution (e.g. for large or high-quality plots):

```julia
A = GeoArray(WorldClim{BioClim}, 5)
plot(A; max_res=3000)
```


Our first example simply loads a file from disk and plots it.

This netcdf file only has one layer, if it has more we could use `GeoStack`
instead.

```@example nc
using GeoData, NCDatasets, Plots 
url = "https://www.unidata.ucar.edu/software/netcdf/examples/tos_O1_2001-2002.nc";
filename = download(url, "tos_O1_2001-2002.nc");
A = GeoArray(filename)
```

Objects with Dimensions other than `X` an `Y` will produce multi-pane plots.
Here we plot every third month in the first year, just using the regular index:

```@example nc
A[Ti(1:3:12)] |> plot
savefig("four_pane.png") 
# output
```

![Global ocean surface temperatures](four_pane.png)

Now plot the ocean temperatures areound the Americas in the first month of 2001.
Notice we are using lat/lon coordinates and date/time instead of regular indexes:

```@example nc
A[Ti(Near(DateTime360Day(2001, 01, 17))), 
  Y(Between(-60.0, 90.0)), 
  X(Between(190.0, 345.0))] |> plot
savefig("tos_americas.png")
```

![Americas regional ocean surface temperature](tos_americas.png)

Now get the mean over the timespan, then save it to disk, and plot it as a
filled contour:

Other plot functions and sliced objects that have only one `X`/`Y`/`Z` dimension
fall back to generic DimensionalData.jl plotting, which will still correctly
label plot axes.


```@exampe nc
using Statistics
# Take the mean
mean_tos = mean(A; dims=Ti)

# Plot a contour plot
contourf(mean_tos; dpi=300, size=(800, 400))
savefig("mean_tos_contour.png")

# Write the mean values to disk
write("mean_tos.nc", mean_tos)
```

![Mean temperatures](mean_tos_contour.png)

Plotting recipes in DimensionalData.jl are the fallback for GeoData.jl when the
object doesn't have 2 `X`/`Y`/`Z` dimensions, or a non-spatial plot command is
used. So (as a random example) we could plot a transect of ocean surface
temperature at 20 degree latitude :

```@example nc
A[Y(Near(20.0)), Ti(1)] |> plot
savefig("tos_transect.png")
```

![Temperatures at lattitude 20-21](tos_transect.png)


GeoData.jl provides a range of other methods that are being added to over time.
Where applicable these methods read and write lazily to and from disk-based
arrays of common raster file types. These methods also work for entire
`GeoStacks` and `GeoSeries` using the same syntax.


## Objects

### GeoArray

Spatial raster data is essentially just an `Array`. But `GeoArray` wrappers
allow treating them as an array that maintains its spatial index, crs and other
metadata through all transformations. This means the can always be plotted and
written to disk after applying most base Julia methods, and most `broadcast`s.

```@docs
AbstractGeoArray
GeoArray
GeoArray(T::Type{<:RasterDataSources.RasterDataSource}, layer)
```

### GeoStack

Spatial data often comes as a bundle of multiple named arrays, as in netcdf.
`GeoStack` can represent this, or multiple files organised in a similar way.

```@docs
AbstractGeoStack
GeoStack
GeoStack(T::Type{<:RasterDataSources.RasterDataSource})
```

### GeoSeries

A series is an meta-array that holds other files/data that is distributed over
some dimension, often time. These files/data can be `GeoArray`s or `GeoStack`s.

```@docs
AbstractGeoSeries
GeoSeries
GeoSeries(T::Type{<:RasterDataSources.RasterDataSource})
```

### Dimensions

GeoData uses `X`, `Y`, and `Z` dimensions from DimensionalData.jl to represent
spatial directions like longitude, latitude and the vertical dimension, and
subset data with them. `Ti` is used for time, and `Band` represent bands. Other
dimensions can have arbitrary names, but will be treated generically. See
[DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl/) for more
details on how they work. 

```@docs
Band
```

### Index modes

These specify properties of the index associated with e.g. the X and Y
dimension. GeoData.jl additional modes to handle dimensions with projections
with `Projected`, and where the projection is mapped to another projection like
`EPSG(4326)` in `Mapped`, to handle e.g. NetCDF dimensions.


```@docs
GeoData.AbstractProjected
Projected
Mapped
```

## Data sources

GeoData.jl uses a number of backends to load raster data. `GeoArray`, `GeoStack`
and `GeoSeries` will detect which backend to use for you, automatically.

#### GRD

R GRD files can be loaded natively, using Julias `MMap` - which means they are
very fast, but are not compressed. They are always 3 dimensional, and have `Y`,
`X` and [`Band`](@ref) dimensions.

#### NetCDF

NetCDF `.nc` files are loaded using
[NCDatasets.jl](https://github.com/Alexander-Barth/NCDatasets.jl). Layers from
files can be loaded as `GeoArray("filename.nc"; key=:layername)`. Without `key`
the first layer is used. `GeoStack("filename.nc")` will use all netcdf variables
in the file that are not dimensions as layers. 

NetCDF layers can have arbitrary dimensions. Known, common dimension names are
converted to `X`, `Y` `Z`, and `Ti`, otherwise `Dim{:layername}` is used. Layers
in the same file may also have different dimensions.

#### GDAL

All files GDAL can access, such as `.tiff` and `.asc` files, can be loaded,
using [ArchGDAL.jl](https://github.com/yeesian/ArchGDAL.jl/issues). These are
generally best loaded as `GeoArray("filename.tif")`, but can be loaded as
`GeoStack("filename.tif"; layersfrom=Band)`, taking layers from the `Band`
dimension, which is also the default.

#### SMAP

The [Soil Moisture Active-Passive](https://smap.jpl.nasa.gov/) dataset provides
global layers of soil moisture, temperature and other related data, in a custom
HDF5 format. Layers are always 2 dimensional, with `Y` and `X` dimensions.

These can be loaded as multi-layered `GeoStack("filename.h5")`. Individual
layers can be loaded as `GeoArray("filename.h5"; key=:layerkey)`, without `key`
the first layer is used.

```@docs
smapseries
```

### RasterDataSources.jl integration

[RasterDataSources.jl](https://github.com/EcoJulia/RasterDataSources.jl)
standardises the download of common raster data sources, with a focus on
datasets used in ecology and the environmental sciences. RasterDataSources.jl is
tightly integrated into GeoData.jl, so that datsets and keywords can be used
directly to download and load data as a `GeoArray`, `GeoStack`, or `GeoSeries`.

```@example
using GeoData, Plots, Dates
A = GeoArray(WorldClim{Climate}, :tavg; month=June)
plot(A)
savefig("worldclim_june_average_temp.png")
```

![WorldClim June verage temperatures](worldclim_june_average_temp.png)

### Writing file formats to disk

Files can be written to disk in all formats other than SMAP HDF5 using
`write("filename.ext", A)`. See the docs for [`write`](@ref). They can (with
some caveats) be written to different formats than they were loaded in,
providing file-type conversion for spatial data.

Some metadata may be lost in formats that store little metadata, or where
metadata conversion has not been completely implemented.


## Exported functions

GeoData.jl is a direct extension of DimensionalData.jl. See [DimensionalData.jl
docs](https://rafaqz.github.io/DimensionalData.jl/stable/) for the majority of
types and functions that can be used in GeoData.jl.

Functions more specific to geospatial data are included in GeoData.jl, and
listed below.

```@docs
aggregate
aggregate!
boolmask
classify 
classify!
convertmode
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
resample
replace_missing
reproject
slice
subset
trim
warp
```

### File operations

These Base and DimensionalData methods have specific GeoData.jl versions:

```@docs
modify
open
read
read!
write
```

## Internals

```@docs
GeoData.FileArray
GeoData.FileStack
GeoData.GeoDiskArray
```
