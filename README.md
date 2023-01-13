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

## Quick start
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
ti = Ti(DateTime(2001):Month(1):DateTime(2001))
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

## Getting the `loopkup` array from dimensions

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
