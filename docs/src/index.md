# GeoData.jl

```@docs
GeoData
```

## Array

Spatial raster data is essentially an array. These wrappers allow treating them
as an array, but maintaining the spacial index, crs and metadata through all
transformations.

```@docs
geoarray
AbstractGeoArray
MemGeoArray
DiskGeoArray
GeoArray
GeoData.OpenGeoArray
```

## Stack

Spatial data often comes as a bundle of multiple named arrays, as in netcdf.
Stacks can represent this, or multiple files organised in a similar way.

```@docs
stack
AbstractGeoStack
MemGeoStack
GeoStack
DiskGeoStack
DiskStack
```

## Series

A series is an meta-array that holds other files/data that is distributed over
some dimension, often time. These files/data can be `geoarray`s or `stack`s.

```@docs
series
AbstractGeoSeries
GeoSeries
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

```@docs
GRDarray
GRDstack
GRDdimMetadata
GRDarrayMetadata
```

## NetCDF

NetCDF files requires NCDatasets.jl to be imported:

```julia
import NCDatasets
```

Single files can be treated as a array or a stack of arrays. 

```@docs
NCDarray
NCDstack
NCDdimMetadata
NCDarrayMetadata
NCDstackMetadata
```

## GDAL

GDAL requires [ArchGDAL.jl](https://github.com/yeesian/ArchGDAL.jl/issues) to be
imported: 

```julia
import ArchGDAL
```

```@docs
GDALarray
GDALstack
GDALdimMetadata
GDALarrayMetadata
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
using [`SMAPseries`](@ref). Methods like `aggregate` can be done over whole
folders of stacks of data with a single command.

```@docs
SMAPdimMetadata
SMAParrayMetadata
SMAPstackMetadata
SMAPstack
SMAPseries
```

## Helper methods

See [DimensionalData.jl docs](https://rafaqz.github.io/DimensionalData.jl/stable/)
for the majority of types and methods that can be used in GeoData.jl. 
GeoData.jl is a direct extension of DimensionalData.jl.

These methods are specific to GeoData.jl:

```@docs
replace_missing
resample
boolmask
missingmask
aggregate
aggregate!
disaggregate
disaggregate!
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
data
```

Not exported:
```@docs
GeoData.filename
GeoData.childkwargs
```

These Base and DimensionalData methods have specific GeoData.jl version:

```@docs
open
write
cat
copy!
modify
```
