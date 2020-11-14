# GeoData.jl

```@docs
GeoData
```

## Dimensions

```@docs
Lat
Lon
Vert
Band
Projected
Converted
```

## Array

```@docs
AbstractGeoArray
MemGeoArray
DiskGeoArray
GeoArray
```

## Stack

```@docs
AbstractGeoStack
MemGeoStack
GeoStack
DiskGeoStack
DiskStack
SMAPstack
```

## Series

```@docs
AbstractGeoSeries
GeoSeries
SMAPseries
```

## Metadata

```@docs
Metadata
DimMetadata
ArrayMetadata
StackMetadata
```

# Sources

## GRD

R GRD files can be loaded natively. The are always 3 dimensional, and have
[`Lat`](@ref), [`Lon`](@ref) and [`Band`](@ref) dimensions.

If ArchGDAL.jl is loaded, they can have [`usercrs`](@ref) and be 

```@docs
GRDarray
GRDstack
GRDdimMetadata
GRDarrayMetadata
```

## NetCDF

NetCDF files required NCDatasets.jl:

```julia
using NCDatasets
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
using ArchGDAL
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
using HDF5
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
write
cat
copy!
replace_missing
boolmask
missingmask
convertmode
reproject
aggregate
aggregate!
disaggregate
disaggregate!
GeoData.alloc_ag
GeoData.alloc_disag
```

Field access:

```@docs
missingval
crs
usercrs
dimcrs
mappedbounds
mappedval
GeoData.filename
```


