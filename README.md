# GeoData

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://rafaqz.github.io/GeoData.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://rafaqz.github.io/GeoData.jl/dev)
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
- Lazy loading: minimisation of RAM requirements for large datasets
- Automation of multi-file/multi-layer tasks with single line commands
- Ecosystem integration: work as much as possible with existing packages
- Ubiquitous DimensionalData.jl dims and selectors for indexing and dimension
  names, but hidden from custom implementations through the AbstractGeoXX interfaces.
- Automatic detection of dimension order, axis range and data orientation and
  order in format-specific implementations.


## Concepts: arrays, stacks, and series

### AbstractGeoArray

`AbstractGeoArray` wraps an array (or location of an array) and metadata 
about its contents. It may be memory (`GeoArray`) or disk-backed (`NCarray`,
`GDAlarray`). They can be indexed as regular Julia arrays or with
DimensionalData.jl dimensions. They will plot as a heatmap in Plots.jl with correct
coordinates and labels, even after slicing with `getindex` or `view`. `getindex`
on a `AbstractGeoArray` will always return a standard `GeoArray`.

### AbstractGeoStack

These are Dict/NamedTuple like structures that may either contain `NamedTuple`
of `AbstractGeoArray`, string paths that will load `AbstractGeoArray`, or a single
path that points to as a multi-layered stack of arrays. 

The primary purpose is that use and syntax is identical for all cases,
abstracting away data source and simplifying access code. `getindex` on any
`AbstractGeoStack` may return a memory backed standard `GeoArray`, or a disk
base AbstractGeoArray. `geoarray[:somelayer] |> plot` plots the layers array,
while `geoarray[:somelayer, Lon(1:100), Band(2)] |> plot` will plot the
subsetted array directly from disk, without loading the whole array. 

`GeoStack` is the generic memory backed type, while `NCstack`, `GDALstack` and
`SMAPstack` are the currently available format-specific implementations.

### AbstractGeoSeries

These are a high-level `AbstractGeoArray` that hold stacks or arrays of paths
they can be loaded from. `GeoSeries` are indexed with dimensions as with a
`AbstractGeoArray`. This is useful when you have multiple files containing
rasters or stacks of rasters spread over dimensions like time and elevation.
As much as possible, implementations should facilitate loading entire
directories and detecting the dimensions from metadata.

This currently allows `series[Time(Near(DateTime(2001,
1))][:temp][Lat(Between(70, 150)), Lon(Between(-20,20))]
|> plot` and in future it will enable syntax like `series[Time(1:2:10), :temp,
Lat(Between(70, 150)), Lon(Between(-20,20))]` to retrieve data from multiple
files efficiently.

`GeoSeries` is the only implementation, as it includes a field indicating its
child type used if loading stacks or arrays from disk.


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
