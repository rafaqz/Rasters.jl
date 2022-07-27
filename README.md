# Rasters

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://rafaqz.github.io/Rasters.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://rafaqz.github.io/Rasters.jl/dev)
[![CI](https://github.com/rafaqz/Rasters.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/rafaqz/Rasters.jl/actions/workflows/ci.yml)
[![Codecov](https://codecov.io/gh/rafaqz/Rasters.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/rafaqz/Rasters.jl)
[![Aqua.jl Quality Assurance](https://img.shields.io/badge/Aquajl-%F0%9F%8C%A2-aqua.svg)](https://github.com/JuliaTesting/Aqua.jl)

Rasters.jl defines common types and methods for reading, writing and
manipulating rasterized spatial data. 

These currently include raster arrays like GeoTIFF and NetCDF, R grd files, 
multi-layered stacks, and multi-file series of arrays and stacks. 

![EarthEnv HabitatHeterogeneity layers trimmed to Australia](https://rafaqz.github.io/Rasters.jl/stable/trim_example_after.png)

_A RasterStack of EarthEnv HabitatHeterogeneity layers, trimmed to Australia and plotted with Plots.jl_

## Data-source abstraction

Rasters provides a standardised interface that allows many source data types to
be used with identical syntax.

- Scripts and packages building on Rasters.jl can treat `Raster`,
  `RasterStack`, and `RasterSeries` as black boxes.
  - The data could hold GeoTiff or NetCDF files, `Array`s in memory or
    `CuArray`s on the GPU - they will all behave in the same way.
  - `RasterStack` can be backed by a Netcdf or HDF5 file, or a `NamedTuple` of
    `Raster` holding `.tif` files, or all `Raster` in memory.
  - Users do not have to deal with the specifics of spatial file types.
- `Projected` lookups with Cylindrical projections can by indexed using other Cylindrical projections
  by setting the `mappedcrs` keyword on construction. You don't need to know the underlying
  projection, the conversion is handled automatically. This means lat/lon
  `EPSG(4326)` can be used seamlessly if you need that.

## Named dimensions and index lookups

Rasters.jl extends
[DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl) so that
spatial data can be indexed using named dimensions like `X`, `Y` and `Ti` (time)
and e.g. spatial coordinates.

Regions and points can be selected with `a..b`, `At` `Near` and `Contains`.

Dimensions can also be used in most `Base` and `Statistics` methods like `mean`
and `reduce` where `dims` arguments are required. Much of the behaviour is
covered in the [DimensionalData
docs](https://rafaqz.github.io/DimensionalData.jl/stable/).

See [the docs](https://rafaqz.github.io/Rasters.jl/stable) for more details and
examples for Rasters.jl.
