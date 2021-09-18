# GeoData.jl

```@docs
GeoData
```

## Examples

The first example loads a file from disk and plots it.

Load GeoData, and NCDatasets, download file and load it to an array. This netcdf
file only has one layer, if it has more we could use `GeoStack` instead.

```@example nc
using GeoData, NCDatasets, Plots 
url = "https://www.unidata.ucar.edu/software/netcdf/examples/tos_O1_2001-2002.nc";
filename = download(url, "tos_O1_2001-2002.nc");
A = GeoArray(filename)
```
Now plot every third month in the first year, just using the regular index:

```@example nc
A[Ti(1:3:12)] |> plot
savefig("four_pane.png") 
# output
```


![Global ocean surface temperatures](four_pane.png)

Now plot Australia in the first month of 2001. Notice we are using lat/lon coordinates 
and date/time instead of regular indexes:

```@example nc
A[Ti(Near(DateTime360Day(2001, 01, 17))), 
  X(Between(0.0, -50.0)), 
  Y(Between(100.0, 160.0))] |> plot
savefig("tos_aus.png")
```

![Australia regional ocean surface temperature](tos_aus.png)

Now get the mean over the timespan, then save it to disk, and plot it :

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

Plotting recipes in DimensionalData.jl are the fallback for GedData.jl when 
the object doesn't have both `X` and `Y` dimensions. So (as a random example) we 
could plot a transect of ocean surface temperature at 20 degree latitude :

```julia
A[Y(Near(20.0)), Ti(1)] |> plot
savefig("tos_transect.png")
```

![Temperatures at lattitude 20-21](tos_transect.png)


GeoData.jl provides a range of other methods that are being added to over time.
Where applicable these methods read and write lazily to and from disk-based
arrays of common raster file types. These methods also work for entire
`GeoStacks` and `GeoSeries` using the same syntax.


# Common Methods

### Changing extent and resolution

DimensionalDat.jl `Dimension`s and `Selctor`s are the primary way to change the
extent of an object. 

|                      |                                                                    |
| :------------------- | :----------------------------------------------------------------- |
| `At(x)`              | get the index exactly matching the passed in value(s)              |
| `Near(x)`            | get the closest index to the passed in value(s)                    |
| `Where(f::Function)` | filter the array axis by a function of the dimension index values. |
| `Between(a, b)`      | get all indices between two values, excluding the high value.      |
| `Contains(x)`        | get indices where the value x falls within an interval             |


Use the `Between` selector to take a `view` of madagascar:

```
A = GeoArray(WorldClim{BioClim}, 5)
madagascar = view(A, X(Between(43.25, 50.48)), Y(Between(-25.61, -12.04))) 
plot(madagascar)
savefig("madagascar_bio5.png")
```

![Temperatures at lattitude 20-21](madagascar_bio5.png)

### Methods that change an objects resolution or extent

|                        |                                                                              |
| :--------------------- | :--------------------------------------------------------------------------- |
| [`aggregate`](@ref)    | aggregate data by the same or different amounts for each axis.               |
| [`disaggregate`](@ref) | similarly disaggregate data.                                                 |
| [`mosaic`](@ref)       | join rasters covering different extents into a single array or file.         |
| [`crop`](@ref)         | shrink objects to specific dimension sizes or the extent of another object.  |
| [`extend`](@ref)       | extend objects to specific dimension sizes or the extent of another object.  |
| [`trim`](@ref)         | trims areas of missing values across arbitrary dimensions and stack layers.  |
| [`resample`](@ref)     | resample data to a different size and projection, or snap to another object. |
| [`warp`](@ref)         | use `gdalwarp` directly on any object, e.g. a multidimensional NetCDF stack. |


### Methods that change an objects values: 

(note that most regular Julia methods work as for any other array)

|                           |                                                                              |
| :------------------------ | :--------------------------------------------------------------------------- |
| [`classify`](@ref)        | classify values into categories.                                             |
| [`mask`](@ref)            | mask and object by a polygon or `GeoArray` along `X/Y`, or other dimensions. |
| [`replace_missing`](@ref) | replace all missing values in an object and update `missingval`.             |


### Methods that change modify data sources:

|                           |                                                                        |
| :------------------------ | :--------------------------------------------------------------------- |
| [`modify`](@ref)          | replace the data in objects. Useful to e.g. move objects to/from a GPU |
| [`read`](@ref)            | read data to memory if it is on disk                                   |
| [`read!`](@ref)           | read data to predefined memory                                         |
| [`open`](@ref)            | open the underlying data for manually reading or writing               |
| [`write`](@ref)           | write objects to file                                                  |


### Altering and summarising arrays and stacks with regular julia methods

Most base methods work as for regular julia `Array`s, from `reverse`, 
rotations like `rotl90`. Base, statistics and linear algebra methods like `mean`
that take a `dims` argument can use the dimension name.

To take the mean over the time dimension:

```julia
mean(A; dims=Ti)
```


`broadcast` works lazily from disk, and is only applied when data is directly
indexed. Adding a dot to any function will use broadcast over a `GeoArray`. 

### Broadcasting

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

`set` can be used to modify the nested properties of dimensions, that are more
difficult to change with `rebuild`. Set works on the principal that dimension
properties can only be in one field, so we generally don't have to specify which
one it is. `set` will also try to update anything affected by a change you make.

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


## Objects

### Array

Spatial raster data is essentially just an `Array`. But `GeoArray` wrappers
allow treating them as an array, while maintaining their spatial index, crs and
other metadata through all transformations. This means the can always be plotted
and written to disk after applying most base Julia methods, and most `broadcast`s.

```@docs
AbstractGeoArray
GeoArray
GeoArray(T::Type{<:RasterDataSources.RasterDataSource}, layer)
```

### Stack

Spatial data often comes as a bundle of multiple named arrays, as in netcdf.
`GeoStack` can represent this, or multiple files organised in a similar way.

```@docs
AbstractGeoStack
GeoStack
GeoStack(T::Type{<:RasterDataSources.RasterDataSource})
```

### Series

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
index into them. `Ti` is used for time, and `Band` represent bands. Other
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

## Sources

GeoData.jl defines a number of wrappers for various file sources. These may
require importing a package before they are available. Increasingly, these
details are not important to the user: `GeoArray`, `GeoStack` and `GeoSeries`
will detect which backend to use for you, automatically.

### GRD

R GRD files can be loaded natively. The are always 3 dimensional, and have
`Y`, `X` and [`Band`](@ref) dimensions.

If ArchGDAL.jl is loaded (to enable reprojection), they can have [`mappedcrs`](@ref).

### NetCDF

NetCDF files are loaded using
[NCDatasets.jl](https://github.com/Alexander-Barth/NCDatasets.jl). Single files
can be treated as GeoArray with a single layer, or a `GeoStack` multiple layers. 

### GDAL

All files GDAL can access can be loaded with GeoData.jl, using
[ArchGDAL.jl](https://github.com/yeesian/ArchGDAL.jl/issues).

### SMAP

The [Soil Moisture Active-Passive](https://smap.jpl.nasa.gov/) dataset provides
global layers of soil moisture, temperature and other related data.

```@docs
smapseries
```


## Plotting

Plots.jl is fully supported. `plot` will plot a heatmap with axes matching
dimension values. If `mappedcrs` is used, converted values will be shown on axes
instead of the underlying `crs` values.

Pixel resolution is limited to allow loading very large files quickly. `max_res` 
specifies the maximum pixel resolution to show on the longest axis of the array.
It can be set manually to change the resolution (e.g. for large or high-quality plots):

```julia
A = GeoArray(WorldClim{BioClim}, 5)
plot(A; max_res=3000)
```

Dimensions other than `X` an `Y` will produce multi-pane plots, while stack
layers will also be plotted in multiple panes, taking the first possible
2-dimensional `X`/`Y`/`Z` slice, along all other dimensions.


## Public interface

GeoData.jl is a direct extension of DimensionalData.jl. See [DimensionalData.jl
docs](https://rafaqz.github.io/DimensionalData.jl/stable/) for the majority of
types and methods that can be used in GeoData.jl.

```@docs
aggregate
aggregate!
boolmask
chunk
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
