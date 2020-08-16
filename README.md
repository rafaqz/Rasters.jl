# GeoData

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://rafaqz.github.io/GeoData.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://rafaqz.github.io/GeoData.jl/dev)
[![Build Status](https://travis-ci.org/rafaqz/GeoData.jl.svg?branch=master)](https://travis-ci.org/rafaqz/GeoData.jl)
[![Codecov](https://codecov.io/gh/rafaqz/GeoData.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/rafaqz/GeoData.jl)

GeoData.jl extends
[DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl) so that
spatial data can be indexed using named dimensions like `Lat` and `Lon`, which
can also be used in most methods like `mean` and `reduce` where dimensions are
required. Much of the behaviour is covered in the [DimensionalData
docs](https://rafaqz.github.io/DimensionalData.jl/stable/).

GeoData.jl also defines common types and methods for accessing and working with
spatial data, such as 2 or multidimensional raster arrays, multi-array "stacks",
and "series" of stacks or arrays that behave similarly or identically across
multiple file-types.

It provides general types `GeoArray`, `GeoStack`, and `GeoSeries`, R `.grd`
files can be loaded natively using `GRDarray` and `GRDstack`. 

GDAL files can be loaded when
[ArchGDAL.jl](https://github.com/yeesian/ArchGDAL.jl) (v0.5 or higher) is
present, with ` GDALarray` and GDALstack. NetCDF similarly can be loaded when
[NCDatasets.jl](https://github.com/Alexander-Barth/NCDatasets.jl) is loaded,
with `NCDarray` and `NCDstack`.

When HDF5.jl is loaded, files from the Soil Moisture Active Passive
([SMAP](https://smap.jpl.nasa.gov/)) dataset can be loaded using `SMAPstack`
or `SMAPseries` to load whole directories. This is both useful for users of
SMAP, and a demonstration of the potential to build standardised interfaces 
for custom spatial dataset formats like those used in SMAP.

Files can be written to disk in all formats using `write`.

Some helper methods for common manipulations are included:
- `aggregate`
- `disaggregate`
- `boolmask`
- `missingmask`
- `replace_missing`

These will be expanded to include interpolation and other tools over time.


**This is a work in progress and the API will break occasionally**

Notably GeoData will shift to relying on
[DiskArrays.jl](https://github.com/meggart/DiskArrays.jl) for wrapping
disk-based data sources when it fully supports GDAL and NetCDF. Currently
this is handled internally. Broadcasting over disk base arrays like `NCDarray`
will be incredibly slow, as chunk-based loading is not implemented.
DiskArrays.jl will solve this, and other problems.

There are no also no guarantees on the accuracy of any of the included methods.
If you are using this in critical applications, please do your own testing,
or add additional tests to the GeoData.jl test suit to verify correctness.

Also note: writing directly to files with `setindex!` is not yet supported. 
You must load to a `GeoArray`, make modifications and `write` a new file. 


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
url = "https://www.unidata.ucar.edu/software/netcdf/examples/tos_O1_2001-2002.nc"
filename = download(url, "tos_O1_2001-2002.nc")
A = NCDarray(filename)
```

Now plot every third month in the first year, just using the regular index:

```julia
using Plots
A[Ti(1:3:12)] |> plot
```

![Global ocean surface temperatures](https://raw.githubusercontent.com/rafaqz/GeoData.jl/media/four_pane_map.png)

Now plot Australia in the first month of 2001.

```julia
A[Ti(Contains(DateTime360Day(2001, 01, 17))), 
  Lat(Between(0, -50)), 
  Lon(Between(100, 160))] |> plot
```

![Australia regional ocean surface temperature](https://raw.githubusercontent.com/rafaqz/GeoData.jl/media/aus.png)

Now get the mean over the timespan, then save it to disk, and plot it :

```julia
using Statistics
mean_tos = mean(A; dims=Ti)
write("mean.ncd, NCDarray, mean_tos))
plot(mean_tos; color=:viridis) 
```

![Mean temperatures](https://raw.githubusercontent.com/rafaqz/GeoData.jl/media/mean.png)

Or a plot transect of ocean surface temperature along the 20 degree latitude line:

```julia
A[Lat(Contains(20)), Ti(1)] |> plot
```

![Temperatures at lattitude 20-21](https://raw.githubusercontent.com/rafaqz/GeoData.jl/media/lat_20.png)


## Works in progress
- Standardised handling and conversion of spatial metadata between data formats
- Handling complex projections: Affine transformation of dimensions to indices.
  AffineMaps will be stored as a wrapper dimension in `dims`.
- Load and write NetCDF projection format
- Integration with Vector/DataFrame spatial types and point/line/polygon data
  types. It should be possible to select polygons of data, and convert between
  linear datasets and array formats.
