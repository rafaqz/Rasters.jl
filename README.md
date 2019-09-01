# GeoData

[![Build Status](https://travis-ci.com/rafaqz/GeoData.jl.svg?branch=master)](https://travis-ci.com/rafaqz/GeoData.jl)
[![Codecov](https://codecov.io/gh/rafaqz/GeoData.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/rafaqz/GeoData.jl)

The core goal of GeoData is to define common types and methods for accessing and 
working with spatial data in Julia, such as 2 or multidimensional raster arrays. 
It provides basic concrete data types, but may also be used as a library to add 
standardised data manipulation and plotting to other geospatial data packages.

GeoData extends
[DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl) so that data
can be indexed using named dimensions, which can also be used in most methods like
`mean` and `reduce` where dimensions are required.

## Design Goals

- Standardisation: data from multiple sources behaves in identical ways, as
  similar to Base methods as possible 
- Easy no-config plotting
- Lazy loading: minimisation of ram requirements for large datasets
- Automation of multi-file/multi-layer tasks with single line commands
- Ecosystem integration: work as much as possible with existing geospatial packages
- Ubiquitous DimensionalData.jl dims and selectors for indexing and dimension
  names, but hidden from custom implementations through the AbstractGeoXX interfaces.
- Automatic detection of dimension order, axis range and data orientation and
  order in implementations


## Concepts: arrays, stacks, and series

### AbstractGeoArray

AbstractGeoArray may be memory (`GeoArray`) or disk backed (`NCarray`,
`GDAlarray`). They can be indexed as regular julia arrays or with
DimensionalData.jl dimensions. They will plot as a heatmap in Plots.jl with correct
coordinates and labels, even after slicing with `getindex` or `view`.

### AbstractGeoStack

These are Dict/NamedTuple like structures that may either contain `NamedTuple`
of `AbstractGeoArray` or paths that will load `AbstractGeoArray`, or a single
path that points to as a multi-layered stack of raster arrays. Use and syntax is
identical for all cases. `geoarray[:somelayer] |> plot` plots the whole array,
while `geoarray[:somelayer, Lon(1:100), Band(2)] |> plot` will plot the
subsetted array directly from disk where possible. 

`GeoStack` is the generic memory backed type, while `NCstack`, `GDALstack` and
`SMAPstack` are the format specific implementations currently available.

### AbstractGeoSeries

These are a high-level `AbstractGeoArray` that hold other stacks or arrays or
point to files that they can be loaded from. `GeoSeries` can be indexed with
dimensions as with a `AbstractGeoArray`. This is useful when you have multiple
files containing rasters or stacks of rasters spread over dimensions like time
and elevation.

This allow `series[Time(Near(DateTime(2001, 1))][:temp][Lat(Between(70, 150)), Lon(Between(-20,20))]
|> plot` and in future will enable syntax like `series[Time(1:2:10), :temp,
Lat(Between(70, 150)), Lon(Between(-20,20))]` to retrieve data from multiple
files efficiently without intermediates.

GeoSeries is the only implementation, which includes a field indicating its
child type used for loading stacks or arrays from disk.



## Works in progress
- Standardised handling of metadata
- Saving data back to disk formats. This should eventually facilitate single
  command conversions between formats, even across whole directories using stack
  and series types.
- Handling complex projections: Affine transformation of dimensions to indices.
  AffineMaps will be stored as a wrapper dimension in `dims`.
- Integration with Vector/DataFrame spatial types and point/line/polygon data
  types. It should be possible to select polygons of data, and convert between
  linear datasets and array formats.
