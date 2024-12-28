# Data sources

Rasters.jl uses a number of backends to load raster data. `Raster`, `RasterStack`
and `RasterSeries` will detect which backend to use for you, automatically.

## GDAL

All files GDAL can access, such as `.tiff` and `.asc` files, can be loaded,
using [ArchGDAL.jl](https://github.com/yeesian/ArchGDAL.jl). These are
generally best loaded as `Raster("filename.tif")`, but can be loaded as
`RasterStack("filename.tif"; layersfrom=Band)`, taking layers from the `Band`
dimension, which is also the default.

## NetCDF

NetCDF `.nc` and some HDF5 `.h5` files cab be loaded using
[NCDatasets.jl](https://github.com/Alexander-Barth/NCDatasets.jl). Layers from
files can be loaded as `Raster("filename.nc"; name=:layername)`. Without `name`
the first layer is used. `RasterStack("filename.nc")` will use all netcdf variables
in the file that are not dimensions as layers. 

NetCDF layers can have arbitrary dimensions. Known, common dimension names are
converted to `X`, `Y` `Z`, and `Ti`, otherwise `Dim{:layername}` is used. Layers
in the same file may also have different dimensions.

## Zarr

Zarr files can be loaded with the [ZarrDatasets.jl](https://github.com/JuliaGeo/ZarrDatasets.jl) 
backend. `Raster(filename; source=Zarrsource())` may be needed where the file type cant be detected 
from the filename. `write` does not yet work for Zarr but will in future.

## GRIB

GRIB files can be loaded with the [ZarrDatasets.jl](https://github.com/JuliaGeo/GRIBDatasets.jl).
`write` is not implemented for GRIB.

## GRD

R GRD files can be loaded natively, using Julias `MMap` - which means they are very fast, but are not compressed. 
They are always 3 dimensional, and have `Y`, `X` and [`Band`](@ref) dimensions.

````@example data_sources
using Rasters
````

````@docs
smapseries
````

## Writing file formats to disk

Files can be written to disk with ArchGDAL.jl and NCDatasets.jl backends using
`write("filename.ext", raster)`. See the docs for [`write`](@ref). 

They can (with some caveats) be written to different formats than they were loaded in as,
providing file-type conversion for spatial data. Some metadata may be lost in formats that 
store little metadata, or where metadata conversion has not been completely implemented.

## RasterDataSources.jl integration

[RasterDataSources.jl](https://github.com/EcoJulia/RasterDataSources.jl)
standardises the download of common raster data sources, with a focus on
datasets used in ecology and the environmental sciences. RasterDataSources.jl is
tightly integrated into Rasters.jl, so that datsets and keywords can be used
directly to download and load data as a `Raster`, `RasterStack`, or `RasterSeries`.

````@example
using Rasters, CairoMakie, Dates
using RasterDataSources
A = Raster(WorldClim{Climate}, :tavg; month=June)
Makie.plot(A)
````

See the docs for [`Raster`](@ref), [`RasterStack`](@ref) and [`RasterSeries`](@ref),
and the docs for `RasterDataSources.getraster` for syntax to specify various
data sources.
