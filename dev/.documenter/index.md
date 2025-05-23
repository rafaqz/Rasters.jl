---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "Rasters.jl"
  text: "Manipulating spatial data"
  tagline: "a powerful package that simplifies the handling of rasterized spatial data in Julia."
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
  - title: 💥💠 Core Functionality
    details: Defines common types and methods for reading, writing, and manipulating rasterized spatial data. <a class="highlight-link">Rasters.jl</a> provides unified data handling through types like <a class="highlight-link">Raster</a>, <a class="highlight-link">RasterStack</a>, and <a class="highlight-link">RasterSeries</a>, offering seamless abstraction regardless of storage backend.
    link: /manual/methods
  - title: ⚙️ Data-Source Abstraction
    details: Rasters provides a standardized interface that enables various data source types to be used with consistent syntax. The data can include <a class="highlight-link">GeoTIFF</a> or <a class="highlight-link">NetCDF</a> files, in-memory Arrays, or <a class="highlight-link">CuArrays</a> on the <a class="highlight-link">GPU</a>—all of which will behave in the same way.
    link: /manual/data_sources
  - title: 🗂️🌐 Data Formats
    details: A <a class="highlight-link">RasterStack</a> can be backed by a <a class="highlight-link">NetCDF</a> or <a class="highlight-link">HDF5</a> file, a NamedTuple of Rasters holding <a class="highlight-link">.tif</a> files, or all Rasters in memory. Users do not need to worry about the specifics of spatial file types.
  - title: 🌍🔍 Effortless Spatial Lookups
    details: <a class="highlight-link">Projected</a> lookups with cylindrical projections can be indexed using other cylindrical projections by setting the <a class="highlight-link">mappedcrs</a> keyword during construction. You don’t need to know the underlying projection, as the conversion is handled automatically. This means that <a class="highlight-link">lat/lon EPSG(4326)</a> can be used seamlessly if needed.
  - title: 🧩⚡ Modular Extensions with Faster Loading
    details: From version 0.8, Rasters.jl uses Julia's package extension system. On Julia 1.9, packages can be placed in <a class="highlight-link">extensions</a> to load only when needed, reducing loading time to under a second. However, each backend or feature <a class="highlight-link">must be loaded manually</a>.
    link: /#Package-extensions
  - title: 🐞📌 Bug reports
    details: Raster data is complex, and there are many opportunities for subtle or obvious bugs to appear. We need bug reports to reduce how often these issues occur over time. However, it's also crucial that issues are easy to reproduce; otherwise, it becomes impractical to fix them. Please follow <a href="#Bugs" class="highlight-link">these guidelines</a>!
    link: /#Bugs
---


## How to Install Rasters.jl? {#How-to-Install-Rasters.jl?}

Since `Rasters.jl` is registered in the Julia General registry, you can simply run the following command in the Julia REPL:

```julia
julia> using Pkg
julia> Pkg.add("Rasters.jl")
# or
julia> ] # ']' should be pressed
pkg> add Rasters
```


If you want to use the latest unreleased version, you can run the following command:

```julia
pkg> add Rasters#main
```


## Package extensions {#Package-extensions}

Before using a specific backend or feature, you must install and load its corresponding package,

::: code-group

```julia [ gdal ]
using Pkg
Pkg.add("ArchGDAL")
```


```julia [ netcdf ]
using Pkg
Pkg.add("NCDatasets")
```


```julia [ grd ]
# built-in
```


```julia [ smap ]
using Pkg
Pkg.add("HDF5")
```


```julia [ grib ]
using Pkg
Pkg.add("GRIBDatasets")
```


```julia [ Makie plots ]
using Pkg
Pkg.add(["GLMakie"])
```


```julia [ Data downloads ]
using Pkg
Pkg.add(["RasterDataSources"]) # which includes Worldclim{Climate}
```


```julia [ Coordinate transforms ]
using Pkg
Pkg.add(["CoordinateTransformations"]) # transformations for gdal rasters
```


:::

and as an example, to use the GDAL backend and download files, you will need to do:

```julia
using Rasters, ArchGDAL, RasterDataSources
```


## 🐞📌  Bugs {#Bugs}

::: info Bugs, errors and making issues for Rasters.jl

Because there are so many raster file types and variations of them, most of the time we need the `exact file` that caused your problem to know how to fix it, and be sure that we have actually fixed it when we are done. So fixing a Rasters.jl bug nearly always involves downloading some file and running some code that breaks with it (if you can trigger the bug without a file, that&#39;s great! but its not always possible).

To make an issue we can fix quickly (or at all) there are three key steps:
1. Include the file in an accessible place on web `without authentication` or any other work on our part, so we can just get it and find your bug. You can put it on a file hosting platform (e.g. google drive, drop box, whatever you use) and share the url.
  
1. Add a [`minimum working example`](https://discourse.julialang.org/t/please-read-make-it-easier-to-help-you/14757) to the issue template that first downloads the file, then runs the function that triggers the bug.
  
1. Paste the `complete stack trace` of the error it produces, right to the bottom, into the issue template. Then we can be sure we reproduced the same problem.
  

Good issues are really appreciated, but they do take just a little extra effort with Rasters.jl because of this need for files.

:::
