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
GDALarray
GrdArray
NCDarray
```

## Stack

```@docs
AbstractGeoStack
MemGeoStack
GeoStack
DiskGeoStack
DiskStack
GDALstack
GrdStack
NCDstack
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
GrdDimMetadata
GrdArrayMetadata
GDALdimMetadata
GDALarrayMetadata
NCDdimMetadata
NCDarrayMetadata
NCDstackMetadata
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
userbounds
userval
GeoData.filename
```


