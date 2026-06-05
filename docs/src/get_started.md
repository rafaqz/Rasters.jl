# Quick start

## Installation

In a script or notebook, install Rasters.jl and other packages used in this tutorial:

````julia
using Pkg
Pkg.add(["Rasters", "Dates", "RasterDataSources", "CairoMakie", "ArchGDAL"])
````

## Creating a Raster

Using Rasters to read GeoTiff or NetCDF files will output something similar to the following toy examples.
This is possible because Rasters.jl extends [DimensionalData.jl](https://rafaqz.github.io/DimensionalData.jl/stable/) so that spatial data can be indexed using named dimensions like `X` and `Y` (coordinates) and `Ti` (time).

````@example first_raster
using Rasters, Dates

lon, lat = X(25:1:30), Y(25:1:30)
ti = Ti(DateTime(2001):Month(1):DateTime(2002))
ras = Raster(rand(lon, lat, ti)) # this generates random numbers with the dimensions given
````

## Getting the lookup array from dimensions

Lookups are ranges assigned to each dimension axis, letting us refer to data within our Raster by its real-world values (like longitude, latitude, or time) instead of integer index positions. 

````@example first_raster
lon = lookup(ras, X) # if X is longitude
lat = lookup(ras, Y) # if Y is latitude
````

!!! info "Dimensions"
      Rasters uses X, Y, and Z dimensions from [`DimensionalData`](https://rafaqz.github.io/DimensionalData.jl/) to represent spatial directions like longitude,
      latitude and the vertical dimension, and subset data with them. Ti is used for time, and Band represents bands.
      Other dimensions can have arbitrary names, but will be treated generically.
      See [`DimensionalData`](https://rafaqz.github.io/DimensionalData.jl/) for more details on how they work.



## Select by index

Select a time slice by its integer index: 

````@example first_raster
ras[Ti(1)]
````
The same slice can also be selected with keyword syntax:

````@example first_raster
ras[Ti=1]
````

A range of indices can also be selected using the syntax `= a:b` or `(a:b)`:

````@example first_raster
ras[Ti(1:10)]
````

## Select by value

Instead of integer positions, we can select using the actual coordinate values with `At`:

````@example first_raster
ras[Ti=At(DateTime(2001))]
````

!!! info "Lookup Arrays"
     These specify properties of the index associated with e.g. the X and Y
     dimension. Rasters.jl defines additional lookup arrays: [`Projected`](@ref) to handle
     dimensions with projections, and [`Mapped`](@ref) where the projection is mapped to
     another projection like `EPSG(4326)`. `Mapped` is largely designed to handle
     NetCDF dimensions, especially with `Explicit` spans.

## Subsetting an object

Regular `getindex` (e.g. `A[1:100, :]`) and `view` work on all objects, just as
with an `Array`. `view` is always lazy, and reads from disk are deferred until
`getindex` is used. 

[`DimensionalData.jl`](https://rafaqz.github.io/DimensionalData.jl/stable/)'s `Dimension`s and `Selector`s are another way to subset an object, using the object's lookups to find values at specific `X/Y` coordinates (or any other dimensions). The available selectors are listed here:

| Selectors |              Description                                                           |
| :--------------------- | :----------------------------------------------------------------- |
| `At(x)`                | get the index exactly matching the passed in value(s).             |
| `Near(x)`              | get the closest index to the passed in value(s).                   |
| `Where(f::Function)`   | filter the array axis by a function of the dimension index values. |
| `a..b`/`Between(a, b)` | get all indices between two values, excluding the high value.      |
| `Contains(x)`          | get indices where the value x falls within an interval.            |

The next example uses some real world data. To download data you will need to specify a folder to put it in. You can do this by assigning the environment variable RASTERDATASOURCES_PATH: 

````julia
ENV["RASTERDATASOURCES_PATH"] = "/home/user/Data/" # your path goes here
````

Use the `..` selector to take a view of Madagascar:

````@example first_raster
import ArchGDAL
using Rasters, RasterDataSources
const RS = Rasters
using CairoMakie
CairoMakie.activate!()

A = Raster(WorldClim{BioClim}, 5)
madagascar = view(A, X(43.25 .. 50.48), Y(-25.61 .. -12.04)) 
# Note the space between .. -12
Makie.plot(madagascar)
````