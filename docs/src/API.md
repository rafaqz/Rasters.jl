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

## Exported functions

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
convertlookup
crop
crs
disaggregate
disaggregate!
extend
extract
inpolygon
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
slice
subset
trim
warp
zonal
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

## Index

```@index
Pages = ["API.md"]
```