# GeoData

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://rafaqz.github.io/GeoData.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://rafaqz.github.io/GeoData.jl/dev)
[![Build Status](https://travis-ci.com/rafaqz/GeoData.jl.svg?branch=master)](https://travis-ci.com/rafaqz/GeoData.jl)
[![Codecov](https://codecov.io/gh/rafaqz/GeoData.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/rafaqz/GeoData.jl)
[![Aqua.jl Quality Assurance](https://img.shields.io/badge/Aquajl-%F0%9F%8C%A2-aqua.svg)](https://github.com/JuliaTesting/Aqua.jl)

GeoData.jl defines common types and methods for reading, writing and
manipulating spatial data, currently raster arrays like GeoTIFF and NetCDF, R
grd files, and multi-array stacks and series of arrays and stack. 


## Lazyness

- Data is loaded lazily wherever possible using DiskArrays.jl. Indexing a
  `GeoStack` by name is always lazy, while `view` of a `GeoArray` is lazy and
  `getindex` will load to memory. `read` can be used on any object to ensure
  that all data is loaded to memory.
- Broadcast over disk-based objects is lazy - it will only run when the array is
  indexed. Always prefer broadcasts to explicit loops - these can be very slow
  with disk-based data.
- Regions and points selected with `Between` and `Contains` select the right
  point or whole interval no matter the order of the index or it's position in
  the cell.

## Data source abstraction

GeoData provides a standardised interface that allows many source data types to
be used with identical syntax.

- Scripts and packages building on GeoData.jl can treat `AbstractGeoArray`,
  `AbstractGeoStack`, and `AbstrackGeoSeries` as black boxes.
  - The data could hold GeoTiff or NetCDF files, `Array`s in memory or
    `CuArray`s on the GPU - they will all behave in the same way.
  - `AbstractGeoStack` can be a Netcdf or HDF5 file, or a `NamedTuple` of
    `GDALarray` holding `.tif` files, or all `GeoArray` in memory.
  - Users do not have to deal with the specifics of spatial file types.
- For `Projected` mode you can index in any projection you want to by setting
  the `mappedcrs` keyword on construction. You don't even need to know the
  underlying projection, the conversion is handled automatically. This means
  lat/lon `EPSG(4326)` can be used across all sources seamlessly if you need
  that.

## Named dimensions and index lookups

GeoData.jl extends
[DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl) so that
spatial data can be indexed using named dimensions like `X`, `Y` and `Ti`
(time), which can also be used in most `Base` and `Statistics` methods like
`mean` and `reduce` where `dims` arguments are required. Much of the behaviour
is covered in the [DimensionalData
docs](https://rafaqz.github.io/DimensionalData.jl/stable/).

## Objects

GeoData.jl provides general types for holding spatial data: `GeoArray`,
`GeoStack`, and `GeoSeries`. 

- `GeoArray` acts like an array, but has attached spatial metadata and can be
  indexed with e.g. it's lattitude/longitude coordinates.
- `GeoStack` is a stack of `GeoArray` layers that share (at least some)
  dimensions. Indexing with `Symbol` names will return a `GeoArray`, indexing
  with other values will return a new `GeoStack` for the selected region, or a
  `NamedTuple` of single values.
- `GeoSeries` is a `GeoArray` of `GeoArray` or `GeoStack`. 
  A `GeoSeries` is most often a time-series, but can be multi-dimensional as well.


## Backends

GeoData.jl relies on ArchGDAL, NCDatasets, HDF5, and MMap for loading
disk-based files. The backend is detected automatically from the file type.

- Netcdf `.nc` files can be loaded with `GeoArray(filename)` or
  `GeoStack(filename)` using 
  [NCDatasets.jl](https://github.com/Alexander-Barth/NCDatasets.jl).
- Anything GDAL can load can be loaded with `GeoArray(filename)`
  [ArchGDAL.jl](https://github.com/yeesian/ArchGDAL.jl). `GeoStack` can be made
  with a `NamedTuple` holding filename `String`s or `GeoArray`s for each stack
  layer.
- `.grd/.gri` files from R can be read and written natively using memory
  mapping. These are usually the most efficient file to work with, but are not
  compressed so have large file sizes.
- Custom HDF5 `.h5` files can also be used. Currently, files from the Soil
  Moisture Active Passive ([SMAP](https://smap.jpl.nasa.gov/)) dataset can be
  loaded with `GeoStack` or `GeoSeries`. This is both useful for users of SMAP,
  and a demonstration of the potential to build standardised interfaces for
  custom spatial formats.

Files can be written to disk in all formats other than HDF5 using `write`, and
can (with some caveats, and loss of metadata) be written to different formats
than they were loaded in, providing file-type conversion for spatial data.
