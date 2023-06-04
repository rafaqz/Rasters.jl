# Rasters

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://rafaqz.github.io/Rasters.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://rafaqz.github.io/Rasters.jl/dev)
[![CI](https://github.com/rafaqz/Rasters.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/rafaqz/Rasters.jl/actions/workflows/ci.yml)
[![Codecov](https://codecov.io/gh/rafaqz/Rasters.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/rafaqz/Rasters.jl)
[![Aqua.jl Quality Assurance](https://img.shields.io/badge/Aquajl-%F0%9F%8C%A2-aqua.svg)](https://github.com/JuliaTesting/Aqua.jl)

[Rasters.jl](https://rafaqz.github.io/Rasters.jl/dev) defines common types and methods for reading, writing and
manipulating rasterized spatial data. 

These currently include raster arrays like GeoTIFF and NetCDF, R grd files, 
multi-layered stacks, and multi-file series of arrays and stacks. 

![EarthEnv HabitatHeterogeneity layers trimmed to Australia](https://rafaqz.github.io/Rasters.jl/stable/trim_example_after.png)

_A RasterStack of EarthEnv HabitatHeterogeneity layers, trimmed to Australia and plotted with Plots.jl_

## Packages extensions and Rasters 0.8 and onwards

On Julia 1.9 we can put additional packages in extensions, so the code only loads when
you load a specific package. Rasters.jl was always intended to work like this,
and its finally possible. This reduced package `using` time from many seconds to well under a second.

But, it means you have to manually load packages you need for each backend or additional
functionality.

For example, to use the GDAL backend, and download files, you now need to do:

```julia
using Rasters, ArchGDAL, RasterDataSources
```

where previously it was just `using Rasters`.

Sources and packages needed:
- `:gdal`: `using ArchGDAL`
- `:netcdf`: `using NCDatasets`
- `:grd`: built-in.
- `:smap`: `using HDF5`
- `:grib`: not yet finished.

Other functionality in extensions:
- Raster data downloads, like `Worldclim{Climate}`: `using RasterDataSources`
- Makie plots: `using Makie`
- Coordinate transformations for gdal rasters: `using CoordinateTransformations`

# Quick start
Install the package by typing:

```julia
]
add Rasters
```

```julia
using Rasters
```

Using `Rasters` to read GeoTiff or NetCDF files will output something similar to the
following toy examples. This is possible because Rasters.jl extends
[DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl) so that
spatial data can be indexed using named dimensions like `X`, `Y` and `Ti` (time)
and e.g. spatial coordinates.

```julia
using Rasters, Dates
lon, lat = X(25:1:30), Y(25:1:30)
ti = Ti(DateTime(2001):Month(1):DateTime(2002))
ras = Raster(rand(lon, lat, ti)) # this generates random numbers with the dimensions given
```
```
6×6×13 Raster{Float64,3} with dimensions: 
  X Sampled{Int64} 25:1:30 ForwardOrdered Regular Points,
  Y Sampled{Int64} 25:1:30 ForwardOrdered Regular Points,
  Ti Sampled{DateTime} DateTime("2001-01-01T00:00:00"):Month(1):DateTime("2002-01-01T00:00:00") ForwardOrdered Regular Points
extent: Extent(X = (25, 30), Y = (25, 30), Ti = (DateTime("2001-01-01T00:00:00"), DateTime("2002-01-01T00:00:00")))
missingval: missing
values: [:, :, 1]
     25         26          27          28         29          30
 25   0.9063     0.427328    0.0320967   0.297023   0.0571002   0.891377
 26   0.443494   0.867547    0.350546    0.150155   0.24565     0.711039
 27   0.745673   0.0991336   0.930332    0.893537   0.805931    0.360583
 28   0.512083   0.125287    0.959434    0.354868   0.337824    0.259563
 29   0.253849   0.692209    0.774092    0.131798   0.823656    0.390013
 30   0.334152   0.136551    0.183555    0.941133   0.450484    0.461862
[and 12 more slices...]
```

## Getting the `lookup` array from dimensions

```julia
lon = lookup(ras, X) # if X is longitude
lat = lookup(ras, Y) # if Y is latitude
```
```
Sampled{Int64} ForwardOrdered Regular Points
wrapping: 25:1:30
```

## Select by index

Selecting a time slice by `index` is done via

```julia
ras[Ti(1)]
```
```
6×6 Raster{Float64,2} with dimensions: 
  X Sampled{Int64} 25:1:30 ForwardOrdered Regular Points,
  Y Sampled{Int64} 25:1:30 ForwardOrdered Regular Points
and reference dimensions: 
  Ti Sampled{DateTime} DateTime("2001-01-01T00:00:00"):Month(1):DateTime("2001-01-01T00:00:00") ForwardOrdered Regular Points
extent: Extent(X = (25, 30), Y = (25, 30))
missingval: missing
values:      25         26          27          28         29          30
 25   0.9063     0.427328    0.0320967   0.297023   0.0571002   0.891377
 26   0.443494   0.867547    0.350546    0.150155   0.24565     0.711039
 27   0.745673   0.0991336   0.930332    0.893537   0.805931    0.360583
 28   0.512083   0.125287    0.959434    0.354868   0.337824    0.259563
 29   0.253849   0.692209    0.774092    0.131798   0.823656    0.390013
 30   0.334152   0.136551    0.183555    0.941133   0.450484    0.461862
```

```julia
ras[Ti=1]
```
```
6×6 Raster{Float64,2} with dimensions: 
  X Sampled{Int64} 25:1:30 ForwardOrdered Regular Points,
  Y Sampled{Int64} 25:1:30 ForwardOrdered Regular Points
and reference dimensions: 
  Ti Sampled{DateTime} DateTime("2001-01-01T00:00:00"):Month(1):DateTime("2001-01-01T00:00:00") ForwardOrdered Regular Points
extent: Extent(X = (25, 30), Y = (25, 30))
missingval: missing
values:      25         26          27          28         29          30
 25   0.9063     0.427328    0.0320967   0.297023   0.0571002   0.891377
 26   0.443494   0.867547    0.350546    0.150155   0.24565     0.711039
 27   0.745673   0.0991336   0.930332    0.893537   0.805931    0.360583
 28   0.512083   0.125287    0.959434    0.354868   0.337824    0.259563
 29   0.253849   0.692209    0.774092    0.131798   0.823656    0.390013
 30   0.334152   0.136551    0.183555    0.941133   0.450484    0.461862
```

or and interval of indices using the syntax `=a:b` or `(a:b)`

```julia
ras[Ti(1:10)]
```
```
6×6×10 Raster{Float64,3} with dimensions: 
  X Sampled{Int64} 25:1:30 ForwardOrdered Regular Points,
  Y Sampled{Int64} 25:1:30 ForwardOrdered Regular Points,
  Ti Sampled{DateTime} DateTime("2001-01-01T00:00:00"):Month(1):DateTime("2001-10-01T00:00:00") ForwardOrdered Regular Points
extent: Extent(X = (25, 30), Y = (25, 30), Ti = (DateTime("2001-01-01T00:00:00"), DateTime("2001-10-01T00:00:00")))
missingval: missing
values: [:, :, 1]
     25         26          27          28         29          30
 25   0.9063     0.427328    0.0320967   0.297023   0.0571002   0.891377
 26   0.443494   0.867547    0.350546    0.150155   0.24565     0.711039
 27   0.745673   0.0991336   0.930332    0.893537   0.805931    0.360583
 28   0.512083   0.125287    0.959434    0.354868   0.337824    0.259563
 29   0.253849   0.692209    0.774092    0.131798   0.823656    0.390013
 30   0.334152   0.136551    0.183555    0.941133   0.450484    0.461862
[and 9 more slices...]
```

## Select by value

```julia
ras[Ti=At(DateTime(2001))]
```
```
6×6 Raster{Float64,2} with dimensions: 
  X Sampled{Int64} 25:1:30 ForwardOrdered Regular Points,
  Y Sampled{Int64} 25:1:30 ForwardOrdered Regular Points
and reference dimensions: 
  Ti Sampled{DateTime} DateTime("2001-01-01T00:00:00"):Month(1):DateTime("2001-01-01T00:00:00") ForwardOrdered Regular Points
extent: Extent(X = (25, 30), Y = (25, 30))
missingval: missing
values:      25         26          27          28         29          30
 25   0.9063     0.427328    0.0320967   0.297023   0.0571002   0.891377
 26   0.443494   0.867547    0.350546    0.150155   0.24565     0.711039
 27   0.745673   0.0991336   0.930332    0.893537   0.805931    0.360583
 28   0.512083   0.125287    0.959434    0.354868   0.337824    0.259563
 29   0.253849   0.692209    0.774092    0.131798   0.823656    0.390013
 30   0.334152   0.136551    0.183555    0.941133   0.450484    0.461862
```

More options are available, like  `Near`, `Contains` and `Where`. For more details go [here](https://rafaqz.github.io/Rasters.jl/dev/#Subsetting-an-object).

Dimensions can also be used in most `Base` and `Statistics` methods like `mean`
and `reduce` where `dims` arguments are required. Much of the behaviour is
covered in the [DimensionalData
docs](https://rafaqz.github.io/DimensionalData.jl/stable/).

See [the docs](https://rafaqz.github.io/Rasters.jl/stable) for more details and
examples for Rasters.jl.

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


# Bugs, errors and making issues for Rasters.jl

Raster data is complicated and there are many places for subtle or not-so-subtle bugs to creep in.

We need bug reports to reduce how often they occur over time. But also, we need issues that are easy to reproduce or it isn't practically possible to fix them.

Because there are so many raster file types and variations of them, most of the time we need the *exact file* that caused your problem to know how to fix it, and be sure that we have actually fixed it when we are done. So fixing a Rasters.jl bug nearly always involves downloading some file and running some code that breaks with it (if you can trigger the bug without a file, thats great! but its not always possible).

To make an issue we can fix quickly (or at all) there are three key steps:

1. Include the file in an accessible place on web *without autentication* or any other work on our part, so we can just get it and find your bug. You can put it on a file hosting platform (e.g. google drive, drop box, whatever you use) and share the url.
2. Add a minimum working example to the issue template that first downloads the file, then runs the function that triggers the bug.
3. Paste the complete stack trace of the error it produces, right to the bottom, into the issue template. Then we can be sure we reproduced the same problem.

Good issues are really appreciated, but they do take just a little extra effort with Rasters.jl because of this need for files.
