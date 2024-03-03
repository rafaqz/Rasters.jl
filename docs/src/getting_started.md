# Quick start

## Install the package by typing:

````julia
] 
add Rasters
````

then do

````julia
using Rasters
````

Using Rasters to read GeoTiff or NetCDF files will output something similar to the following toy examples.
This is possible because Rasters.jl extends DimensionalData.jl so that spatial data can be indexed using
named dimensions like `X`, `Y` and `Ti` (time) and e.g. spatial coordinates.

````@example first_raster
using Rasters, Dates

lon, lat = X(25:1:30), Y(25:1:30)
ti = Ti(DateTime(2001):Month(1):DateTime(2002))
ras = Raster(rand(lon, lat, ti)) # this generates random numbers with the dimensions given
````

## Getting the lookup array from dimensions

````@example first_raster
lon = lookup(ras, X) # if X is longitude
lat = lookup(ras, Y) # if Y is latitude
````

## Select by index
Selecting a time slice by index is done via
````@example first_raster
ras[Ti(1)]
````
also

````@example first_raster
ras[Ti=1]
````

or and interval of indices using the syntax =a:b or (a:b)
````@example first_raster
ras[Ti(1:10)]
````
## Select by value

````@example first_raster
ras[Ti=At(DateTime(2001))]
````

More options are available, like `Near`, `Contains` and `Where`.

!!! info "Dimensions"
      Rasters uses X, Y, and Z dimensions from [`DimensionalData`](https://rafaqz.github.io/DimensionalData.jl/) to represent spatial directions like longitude,
      latitude and the vertical dimension, and subset data with them. Ti is used for time, and Band represent bands.
      Other dimensions can have arbitrary names, but will be treated generically.
      See [`DimensionalData`](https://rafaqz.github.io/DimensionalData.jl/) for more details on how they work.


!!! info "Lookup Arrays"
     These specify properties of the index associated with e.g. the X and Y
     dimension. Rasters.jl defines additional lookup arrays: [`Projected`](@ref) to handle
     dimensions with projections, and [`Mapped`](@ref) where the projection is mapped to
     another projection like `EPSG(4326)`. `Mapped` is largely designed to handle
     NetCDF dimensions, especially with `Explicit` spans.


## Subsetting an object

Regular `getindex` (e.g. `A[1:100, :]`) and `view` work on all objects just as
with an `Array`. `view` is always lazy, and reads from disk are deferred until
`getindex` is used. `DimensionalData.jl` `Dimension`s and `Selector`s are the other
way to subset an object, making use of the objects index to find values at 
e.g. certain `X/Y` coordinates. The available selectors are listed here:

| <div style="width:120px">Selectors</div> |              Description                                                           |
| :--------------------- | :----------------------------------------------------------------- |
| `At(x)`                | get the index exactly matching the passed in value(s).             |
| `Near(x)`              | get the closest index to the passed in value(s).                   |
| `Where(f::Function)`   | filter the array axis by a function of the dimension index values. |
| `a..b`/`Between(a, b)` | get all indices between two values, excluding the high value.      |
| `Contains(x)`          | get indices where the value x falls within an interval.            |

!!! info
      - Use the `..` selector to take a view of madagascar:

````@example first_raster
using Rasters, RasterDataSources
const RS = Rasters
using CairoMakie
CairoMakie.activate!()

A = Raster(WorldClim{BioClim}, 5)
madagascar = view(A, X(43.25 .. 50.48), Y(-25.61 .. -12.04)) 
# Note the space between .. -12
Makie.plot(madagascar)
````