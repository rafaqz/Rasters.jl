# GeoData

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://rafaqz.github.io/GeoData.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://rafaqz.github.io/GeoData.jl/dev)
[![Build Status](https://travis-ci.org/rafaqz/GeoData.jl.svg?branch=master)](https://travis-ci.org/rafaqz/GeoData.jl)
[![Codecov](https://codecov.io/gh/rafaqz/GeoData.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/rafaqz/GeoData.jl)

The core goal of GeoData is to define common types and methods for accessing and
working with spatial data in Julia, such as 2 or multidimensional raster arrays.
It provides basic concrete data types, but may also be used as a library to add
standardised data manipulation and plotting to other geospatial data packages.

GeoData.jl extends
[DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl) so that data
can be indexed using named dimensions, which can also be used in most methods
like `mean` and `reduce` where dimensions are required. Most behaviour is
covered in the [DimensionalData docs](https://rafaqz.github.io/DimensionalData.jl/stable/).

## Design Goals

- Standardisation: data from multiple sources have similar or identical behaviour.
- Easy, no-config plotting
- Lazy loading: minimisation of RAM requirements for large datasets
- Automation of multi-file/multi-layer tasks with single line commands
- Ecosystem integration: work as much as possible with existing packages
- Ubiquitous DimensionalData.jl dims and selectors for indexing and dimension
  names, but hidden from custom implementations through the AbstractGeoX interfaces.
- Automatic detection of dimension order, axis range and data orientation and
  order in format-specific implementations.

## Works in progress
- Standardised handling of metadata
- Handling complex projections: Affine transformation of dimensions to indices.
  AffineMaps will be stored as a wrapper dimension in `dims`.
- Integration with Vector/DataFrame spatial types and point/line/polygon data
  types. It should be possible to select polygons of data, and convert between
  linear datasets and array formats.
