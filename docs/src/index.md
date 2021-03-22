# GeoData.jl

```@docs
GeoData
```

## Dimensions

```@docs
Band
```

## Index modes

```@docs
GeoData.AbstractProjected
Projected
Mapped
```

## Array

```@docs
geoarray
AbstractGeoArray
MemGeoArray
DiskGeoArray
GeoArray
GeoData.OpenGeoArray
```

## Stack

```@docs
stack
AbstractGeoStack
MemGeoStack
GeoStack
DiskGeoStack
DiskStack
SMAPstack
```

## Series

```@docs
series
AbstractGeoSeries
GeoSeries
SMAPseries
```

# Sources

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

NetCDF files required NCDatasets.jl:

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
available: 

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

It uses a custom format of HDF5 files, so required HDF5.jl to be available:

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
