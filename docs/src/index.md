# GeoData.jl

```@docs
GeoData
```

## Array

Spatial raster data is essentially an array. These wrappers allow treating them
as an array, but maintaining the spacial index, crs and metadata through all
transformations.

```@docs
AbstractGeoArray
GeoArray
GeoArray(T::Type{<:RasterDataSources.RasterDataSource})
```

## Stack

Spatial data often comes as a bundle of multiple named arrays, as in netcdf.
Stacks can represent this, or multiple files organised in a similar way.

```@docs
AbstractGeoStack
GeoStack
GeoStack(T::Type{<:RasterDataSources.RasterDataSource})
```

## Series

A series is an meta-array that holds other files/data that is distributed over
some dimension, often time. These files/data can be `geoarray`s or `stack`s.

```@docs
AbstractGeoSeries
GeoSeries
GeoSeries(T::Type{<:RasterDataSources.RasterDataSource})
```

## Dimensions

GeoData uses `X`, `Y`, and `Z` dimensions from DimensionalData.jl to represent
spatial directions like longitude, latitude and the vertical dimension, and
index into them. See
[DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl/) for more
details on how they work. GeoData.jl defines a `Band` dimension to represent
bands.

```@docs
Band
```

## Index modes

```@docs
GeoData.AbstractProjected
Projected
Mapped
```

# Sources

GeoData.jl defines a number of wrappers for various file sources. These may
require importing a package before they are available. Increasingly, these
details are not important to the user: `geoarray`, `stack` and `series` will
detect which backend to use for you, automatically.

## GRD

R GRD files can be loaded natively. The are always 3 dimensional, and have
`Y`, `X` and [`Band`](@ref) dimensions.

If ArchGDAL.jl is loaded (to enable reprojection), they can have [`mappedcrs`](@ref).

## NetCDF

NetCDF files requires NCDatasets.jl to be imported:

```julia
import NCDatasets
```

Single files can be treated as a array or a stack of arrays. 


## GDAL

GDAL requires [ArchGDAL.jl](https://github.com/yeesian/ArchGDAL.jl/issues) to be
imported: 

```julia
import ArchGDAL
```

## SMAP

The [Soil Moisture Active-Passive](https://smap.jpl.nasa.gov/) dataset provides
global layers of soil moisture, temperature and other related data.

It uses a custom format of HDF5 files, so requires HDF5.jl to be imported:

```julia
import HDF5
```

Files must be downloaded manually due to authentication restrictions. As the
datasets are know files in standardised formats, whole folders can be loaded
using [`smapseries`](@ref). Methods like `aggregate` can be done over whole
folders of stacks of data with a single command.

```@docs
smapseries
```

## Plotting

Plots.jl is fully supported. `plot` will plot a heatmap with axes matching
dimension values. If `mappedcrs` is used converted values will be shown on 
axes instead of the underlying `crs` values.

Pixel resolution is limited to allow loading very large files. `max_res` 
specifies the maximum pixel resolution to show on the longest axis of the array.
It can be set manually to change the resolution (e.g. for large or high-quality plots):

```julia
plot(A; max_res=3000)
```

Dimensions other than `X` an `Y` will produce multi-pane plots.

![Global ocean surface temperatures](https://raw.githubusercontent.com/rafaqz/GeoData.jl/media/four_pane_map.png)


## Helper methods

See [DimensionalData.jl docs](https://rafaqz.github.io/DimensionalData.jl/stable/)
for the majority of types and methods that can be used in GeoData.jl. 
GeoData.jl is a direct extension of DimensionalData.jl.

These methods are specific to GeoData.jl:

```@docs
replace_missing
boolmask
missingmask
aggregate
aggregate!
disaggregate
disaggregate!
resample
warp
crop
trim
extend
slice
chunk
points
convertmode
reproject
```

Field access:

```@docs
missingval
crs
mappedcrs
mappedbounds
mappedindex
```

These Base and DimensionalData methods have specific GeoData.jl versions:

```@docs
modify
read
read!
open
write
```

## Internals

```@docs
GeoData.FileArray
GeoData.FileStack
GeoData.GeoDiskArray
```
