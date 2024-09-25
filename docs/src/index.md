```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "Rasters.jl"
  tagline: "Manipulating rasterized spatial data"
  image:
    src: /logo.png
    alt: Rasters
  actions:
    - theme: brand
      text: Get Started
      link: /get_started
    - theme: alt
      text: View on Github
      link: https://github.com/rafaqz/Rasters.jl
    - theme: alt
      text: API Reference
      link: /api
      

features:
  - title: Rasters.jl
    details: Defines common types and methods for reading, writing and manipulating rasterized spatial data.
    link: /markdown-examples
  - title: Data Formats
    details: These currently include raster arrays like <strong>GeoTIF</strong> and <strong>NetCDF</strong>, <strong>R grd</strong> files, multi-layered stacks, and multi-file series of arrays and stacks.
---
```


!!! info "Data-source abstraction"

    Rasters provides a standardised interface that allows many source data types to
    be used with identical syntax.

    - Scripts and packages building on Rasters.jl can treat `Raster`, `RasterStack`, and `RasterSeries` as black boxes.
    - The data could hold GeoTiff or NetCDF files, `Array`s in memory or `CuArray`s on the GPU - they will all behave in the same way.
    - `RasterStack` can be backed by a Netcdf or HDF5 file, or a `NamedTuple` of `Raster` holding `.tif` files, or all `Raster` in memory.
    - Users do not have to deal with the specifics of spatial file types.
    - `Projected` lookups with Cylindrical projections can by indexed using other Cylindrical projections by setting the `mappedcrs` keyword on construction. You don't need to know the underlying projection, the conversion is handled automatically. This means lat/lon `EPSG(4326)` can be used seamlessly if you need that.

## Installation

````julia
] add Rasters
````

## Packages extensions

!!! tip "Packages extensions and Rasters 0.8 and onwards"

    On Julia 1.9 we can put additional packages in extensions, so the code only loads when
    you load a specific package. Rasters.jl was always intended to work like this,
    and its finally possible. This reduced package `using` time from many seconds to well under a second.

    But, it means you have to manually load packages you need for each backend or additional
    functionality.

For example, to use the GDAL backend, and download files, you now need to do:

```julia
using Rasters, ArchGDAL, RasterDataSources
```
where previously it was just using Rasters.

Sources and packages needed:

- :gdal: `using ArchGDAL`
- :netcdf: `using NCDatasets`
- :grd: built-in.
- :smap: `using HDF5`
- :grib: not yet finished.

Other functionality in extensions:

- Raster data downloads, like Worldclim{Climate}: `using RasterDataSources`
- Makie plots: `using Makie`
- Coordinate transformations for gdal rasters: `using CoordinateTransformations`

## Bugs and errors

::: warning Bugs, errors and making issues for Rasters.jl

Raster data is complicated and there are many places for subtle or not-so-subtle bugs to creep in.

We need bug reports to reduce how often they occur over time. But also, we need issues that are easy to reproduce or it isn't practically possible to fix them.

Because there are so many raster file types and variations of them, most of the time we need the *exact file* that caused your problem to know how to fix it, and be sure that we have actually fixed it when we are done. So fixing a Rasters.jl bug nearly always involves downloading some file and running some code that breaks with it (if you can trigger the bug without a file, that's great! but its not always possible).

To make an issue we can fix quickly (or at all) there are three key steps:

1. Include the file in an accessible place on web *without authentication* or any other work on our part, so we can just get it and find your bug. You can put it on a file hosting platform (e.g. google drive, drop box, whatever you use) and share the url.
2. Add a minimum working example to the issue template that first downloads the file, then runs the function that triggers the bug.
3. Paste the complete stack trace of the error it produces, right to the bottom, into the issue template. Then we can be sure we reproduced the same problem.

Good issues are really appreciated, but they do take just a little extra effort with Rasters.jl because of this need for files.

:::
