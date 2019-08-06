# GeoData

[![Build Status](https://travis-ci.com/rafaqz/GeoData.jl.svg?branch=master)](https://travis-ci.com/rafaqz/GeoData.jl)
[![Codecov](https://codecov.io/gh/rafaqz/GeoData.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/rafaqz/GeoData.jl)

The core goal of GeoData is to define common types and methods for accessing and 
working with spatial data in Julia, such as 2 or multidimensional raster arrays. 
It provides basic concrete data types, but may also be used as a library to add 
standardised data manipulation and plotting to other geospatial data packages.

## Features

Leveraging DimentionalData.jl:

- Spatial data remains attached to arrays after subsetting with getindex or view, and is updated
  to match the subset where necessary.
- Dimensions can be handled in any order: Lat, Lon, Vert and Time. Custom dims can be added.
- Common plotting recipes work for all AbstractGeoArrays


## Works in progress
- The ability to convert data between different types automatically, for example 
  to extract a matrix from a NetCDF file that keeps information about its 
  coordinates and projection, and plots correctly. This requires packages depending 
  on DimensionalData (at minimum) or GeoData.
- Affine transformation of dimensions to indices. AffineMaps can be stored as a
  wrapper dimension, eg. `AffineDims(Lat(), Lon(), afmap)` so that Lon(x), Lat(y) 
  indexing can be used for directly selecting indices, and affine
  transformations can be used for `select()`.
- Integration with point/line/plygon data types. It should be possible to select
  polygons of data, and convert between linear dataset and array formats.
