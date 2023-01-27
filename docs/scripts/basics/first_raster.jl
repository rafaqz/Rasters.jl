# Rasters.jl defines common types and methods for reading, writing and manipulating rasterized spatial data.
# These currently include raster arrays like `GeoTIFF` and `NetCDF`, `R grd` files, multi-layered stacks, and multi-file series of arrays and stacks. 

# ## Subsetting an object

# Regular `getindex` (e.g. `A[1:100, :]`) and `view` work on all objects just as
# with an `Array`. `view` is always lazy, and reads from disk are deferred until
# `getindex` is used. `DimensionalData.jl` `Dimension`s and `Selector`s are the other
# way to subset an object, making use of the objects index to find values at 
# e.g. certain X/Y coordinates. The available selectors are listed here:

# |                        |                                                                    |
# | :--------------------- | :----------------------------------------------------------------- |
# | `At(x)`                | get the index exactly matching the passed in value(s).             |
# | `Near(x)`              | get the closest index to the passed in value(s).                   |
# | `Where(f::Function)`   | filter the array axis by a function of the dimension index values. |
# | `a..b`/`Between(a, b)` | get all indices between two values, excluding the high value.      |
# | `Contains(x)`          | get indices where the value x falls within an interval.            |

# !!! info
#      - Use the `..` selector to take a view of madagascar:


using Rasters
const RS = Rasters
using CairoMakie
CairoMakie.activate!()

A = Raster(WorldClim{BioClim}, 5)
madagascar = view(A, X(43.25 .. 50.48), Y(-25.61 .. -12.04)) # Note the space between .. -12
#plot(madagascar)
#lookup(madagascar, RS.Band(1))
mad = replace_missing(madagascar[RS.Band(1)], missingval=NaN)
x, y = lookup(mad, X), lookup(mad, Y)

heatmap(x, y, mad.data;
    axis = (; aspect = DataAspect()),
    figure = (; resolution=(600,400))
    )
#current_figure()