# GeoData

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://rafaqz.github.io/GeoData.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://rafaqz.github.io/GeoData.jl/dev)
[![Build Status](https://travis-ci.org/rafaqz/GeoData.jl.svg?branch=master)](https://travis-ci.org/rafaqz/GeoData.jl)
[![Codecov](https://codecov.io/gh/rafaqz/GeoData.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/rafaqz/GeoData.jl)

GeoData.jl defines common types and methods for accessing and
working with spatial data in Julia, such as 2 or multidimensional raster arrays.
It provides general types `GeoArray`, `GeoStack`, and `GeoSeries`, and
source specific types for loading GDAL, NetCDF and other file types, 
available when packages like ArchGDAL.jl or NCDatasets.jl are loaded.

GeoData.jl is useful both as a scripting tool, and as a library of 
standardised data manipulation for use in other geospatial data and
modelling packages.

GeoData.jl extends
[DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl) so that data
can be indexed using named dimensions, which can also be used in most methods
like `mean` and `reduce` where dimensions are required. Most behaviour is
covered in the [DimensionalData docs](https://rafaqz.github.io/DimensionalData.jl/stable/).

## Goals

- Standardisation: data from multiple sources has similar or identical syntax
  and behaviour.
- Easy, no-config plotting
- Lazy loading: minimisation of memory requirements for large datasets
- Accuracy: `Selector`s should select exact regions, and handle points both 
  and intervals. 
- Multi-layer, multi-file objects. `GeoStack` and `GeoSeries` facilitate
  simple operations over large datasets, with detail abstracted away from
  users and other packages.

## Examples

We'll load a file from disk, and do some manipulations and plotting.

Load GeoData, and NCDatasets, download file and load it to 
an array. This netcdf file only has one layer, if it has more we 
could use `NCDstack` instead.

```julia
using GeoData, NCDatasets
filename = download("https://www.unidata.ucar.edu/software/netcdf/examples/tos_O1_2001-2002.nc", "tos_O1_2001-2002.nc")
A = NCDarray(filename)
```

Now plot every third month in the first year, just using the regular index:

```julia
using Plots
A[Ti(1:3:12)] |> plot
```

Now plot Australia in the first month of 2001.

```julia
A[Ti(Contains(DateTime360Day(2001, 01, 17))), Lat(Between(0, -50)), Lon(Between(100, 160))] |> plot
```

Now get the mean over the timespan, then save it to disk, and plot it :

```julia
using Statistics
mean_tos = mean(A; dims=Ti)
write("mean.ncd, NCDarray, mean_tos))
plot(mean_tos; color=:viridis) 
```

Or a transect of ocean surface temperature along the 20 degree latitude line:

```julia
A[Lat(Contains(20)), Ti(1)] |> plot
```
pyplot()
filename = download("https://www.unidata.ucar.edu/software/netcdf/examples/tos_O1_2001-2002.nc", "tos_O1_2001-2002.nc")
A = NCDarray(filename)



write("mean.netcdf", NCDarray, mean(A; dims=Ti))

A[Lat(Contains(20)), Ti(1)] |> plot


## Works in progress
- Standardised handling of metadata
- Handling complex projections: Affine transformation of dimensions to indices.
  AffineMaps will be stored as a wrapper dimension in `dims`.
- Integration with Vector/DataFrame spatial types and point/line/polygon data
  types. It should be possible to select polygons of data, and convert between
  linear datasets and array formats.
