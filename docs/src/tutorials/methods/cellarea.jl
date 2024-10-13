#=
# `cellarea` tutorial

```@meta
CollapsedDocStrings=true
```

```@docs; canonical=false
cellarea
```

`cellarea` computes the area of each cell in a raster.  
This is useful for a number of reasons - if you have a variable like 
population per cell, or elevation ([spatially extensive variables](https://r-spatial.org/book/05-Attributes.html#sec-extensiveintensive)),
you'll want to account for the fact that different cells have different areas.

`cellarea` returns a Raster with the same x and y dimensions as the input, 
where each value in the raster encodes the area of the cell (in meters by default).

You can specify whether you want to compute the area in the plane of your projection 
(`Planar()`), on a sphere of some radius (`Spherical(; radius=...)`).
<!-- or on an ellipsoid 
(`Geodetic()`), using the first argument.-->

Let's construct a raster and see what this looks like!  We'll keep it in memory.

The spherical <!-- and geodetic --> method relies on the [Proj.jl](https://github.com/JuliaGeo/Proj.jl) package to perform coordinate transformation, so that has to be loaded explicitly.
=#

using Rasters
import Proj # to activate the spherical `cellarea` method
import ArchGDAL # purely for data loading

# To construct a raster, we'll need to specify the x and y dimensions.  These are called "lookups" in Rasters.jl.
using Rasters.Lookups
# We can now construct the x and y lookups.  Here we'll use a start-at-one, step-by-five grid.
# Note that we're specifying that the "sampling", i.e., what the coordinates actually mean, 
# is `Intervals(Start())`, meaning that the coordinates are the starting point of each interval.
#
# This is in contrast to `Points()` sampling, where each index in the raster represents the value at a sampling point;
# here, each index represents a grid cell, which is defined by the coordinate being at the start.
x = X(1:5:30; sampling = Intervals(Start()), crs = EPSG(4326))
y = Y(50:5:80; sampling = Intervals(Start()), crs = EPSG(4326))
# I have chosen the y-range here specifically so we can show the difference between spherical and planar `cellarea`.
ras = Raster(ones(x, y); crs = EPSG(4326))
# We can just call `cellarea` on this raster, which returns cell areas in meters, on Earth, assuming it's a sphere:
cellarea(ras)
# and if we plot it, you can see the difference in cell area as we go from the equator to the poles:
using Makie, CairoMakie
heatmap(cellarea(ras); axis = (; aspect = DataAspect()))
# We can also try this using the planar method, which simply computes the area of the rectangle using `area = x_side_length * y_side_length`:
cellarea(Planar(), ras)
# Note that this is of course wildly inaccurate for a geographic dataset - but if you're working in a projected coordinate system, like polar stereographic or Mercator, this can be very useful (and a _lot_ faster)!

#=
## Usage example
Let's get the rainfall over Chile, and compute the average rainfall per meter squared across the country for the month of June.

We'll get the precipitation data across the globe from [WorldClim](https://www.worldclim.org/data/index.html), via [RasterDataSources.jl](https://github.com/EcoJulia/RasterDataSources.jl), and use the `month` keyword argument to get the June data.

Then, we can get the geometry of Chile from [NaturalEarth.jl](https://github.com/JuliaGeo/NaturalEarth.jl), and use `Rasters.mask` to get the data just for Chile.
=#

using RasterDataSources, NaturalEarth

precip = Raster(WorldClim{Climate}, :prec; month = 6)
#
all_countries = naturalearth("admin_0_countries", 10)
chile = all_countries.geometry[findfirst(==("Chile"), all_countries.NAME)]
# Let's plot the precipitation on the world map, and highlight Chile:
f, a, p = heatmap(precip; colorrange = Makie.zscale(replace_missing(precip, NaN)), axis = (; aspect = DataAspect()))
p2 = poly!(a, chile; color = (:red, 0.3), strokecolor = :red, strokewidth = 0.5)
f
# You can see Chile highlighted in red, in the bottom left quadrant.
#
# First, let's make sure that we only have the data that we care about, and crop and mask the raster so it only has values in Chile.
# We can crop by the geometry, which really just generates a view into the raster that is bounded by the geometry's bounding box.
cropped_precip = crop(precip; to = chile)
# Now, we mask the data such that any data outside the geometry is set to `missing`.
masked_precip = mask(cropped_precip; with = chile)
heatmap(masked_precip)
# This is a lot of missing data, but that's mainly because the Chile geometry we have encompasses the Easter Islands as well, in the middle of the Pacific.

# Now, let's compute the average precipitation per square meter across Chile.
# First, we need to get the area of each cell in square meters.  We'll use the spherical method, since we're working with a geographic coordinate system.  This is the default.
areas = cellarea(masked_precip)
masked_areas = mask(areas; with = chile)
heatmap(masked_areas; axis = (; title = "Cell area in square meters"))
# Now we can compute the average precipitation per square meter.  First, we compute total precipitation per grid cell:
precip_per_area = masked_precip .* masked_areas
# We can sum this to get the total precipitation per square meter across Chile:
total_precip = sum(skipmissing(precip_per_area))
# We can also sum the areas to get the total area of Chile (in this raster, at least).
total_area = sum(skipmissing(masked_areas))
# And we can convert that to an average by dividing by the total area:
avg_precip = total_precip / total_area
# According to the internet, Chile gets about 100mm of rain per square meter in June, so our statistic seems pretty close.
#
# Let's see what happens if we don't account for cell areas:
bad_total_precip = sum(skipmissing(masked_precip))
bad_avg_precip = bad_total_precip / length(collect(skipmissing(masked_precip)))
# This is misestimated!  This is why it's important to account for cell areas when computing averages.

#=
!!! note
    If you made it this far, congratulations!  

    It's interesting to note that we've replicated the workflow of `zonal` here.  
    `zonal` is a more general function that can be used to compute any function over geometries, 
    and it has multithreading built in.  
    
    But fundamentally, this is all that `zonal` is doing under the hood - 
    masking and cropping the raster to the geometry, and then computing the statistic.
=#