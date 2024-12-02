# Computing spatial means

```@meta
CollapsedDocStrings=true
```

It's very common to want to compute the mean of some value over some area of a raster.  The initial approach is to simply average the values, but this will give you the arithmetic mean, not the spatial mean.

The reason for this is that raster cells do not always have the same area, especially over a large region of the Earth where its curvature comes into play.

To compute the spatial mean, you need to weight the values by the area of each cell.  You can do this by multiplying the values by the cell area, then summing the values, and dividing that number by the total area.  That was the motivation for this example.

Let's get the rainfall over Chile, and compute the average rainfall across the country for the month of June.

## Acquiring the data

We'll get the precipitation data across the globe from [WorldClim](https://www.worldclim.org/data/index.html), via [RasterDataSources.jl](https://github.com/EcoJulia/RasterDataSources.jl), and use the `month` keyword argument to get the June data.

Then, we can get the geometry of Chile from [NaturalEarth.jl](https://github.com/JuliaGeo/NaturalEarth.jl), and use `Rasters.mask` to get the data just for Chile.

````@example cellarea
using Rasters
import Proj # to activate the spherical `cellarea` method

using ArchGDAL, RasterDataSources, NaturalEarth # purely for data loading

using CairoMakie # for plotting

precip = Raster(WorldClim{Climate}, :prec; month = 6)
````

````@example cellarea
all_countries = naturalearth("admin_0_countries", 10)
chile = all_countries.geometry[findfirst(==("Chile"), all_countries.NAME)]
````

Let's plot the precipitation on the world map, and highlight Chile:

````@example cellarea
f, a, p = heatmap(precip; colorrange = Makie.zscale(replace_missing(precip, NaN)), axis = (; aspect = DataAspect()))
p2 = poly!(a, chile; color = (:red, 0.3), strokecolor = :red, strokewidth = 0.5)
f
````

You can see Chile highlighted in red, in the bottom left quadrant.

## Processing the data

First, let's make sure that we only have the data that we care about, and crop and mask the raster so it only has values in Chile.
We can crop by the geometry, which really just generates a view into the raster that is bounded by the geometry's bounding box.

````@example cellarea
cropped_precip = crop(precip; to = chile)
````

Now, we mask the data such that any data outside the geometry is set to `missing`.

````@example cellarea
masked_precip = mask(cropped_precip; with = chile)
heatmap(masked_precip)
````

This is a lot of missing data, but that's mainly because the Chile geometry we have encompasses the Easter Islands as well, in the middle of the Pacific.


```@docs; canonical=false
cellarea
```

`cellarea` computes the area of each cell in a raster.
This is useful for a number of reasons - if you have a variable like
population per cell, or elevation ([spatially extensive variables](https://r-spatial.org/book/05-Attributes.html#sec-extensiveintensive)),
you'll want to account for the fact that different cells have different areas.

You can specify whether you want to compute the area in the plane of your projection
(`Planar()`), or on a sphere of some radius (`Spherical(; radius=...)`).

Now, let's compute the average precipitation per square meter across Chile.
First, we need to get the area of each cell in square meters.  We'll use the spherical method, since we're working with a geographic coordinate system.  This is the default.

````@example cellarea
areas = cellarea(masked_precip)
masked_areas = mask(areas; with = chile)
heatmap(masked_areas; axis = (; title = "Cell area in square meters"))
````

You can see here that cells are largest towards the equator, and smallest away from it.  This means that cells away from the equator should have a smaller contribution to the average than cells nearer the equator.

## Computing the spatial mean

Now we can compute the average precipitation per square meter.  First, we compute total precipitation per grid cell:

````@example cellarea
precip_per_area = masked_precip .* masked_areas
````

We can sum this to get the total precipitation per square meter across Chile:

````@example cellarea
total_precip = sum(skipmissing(precip_per_area))
````

We can also sum the areas to get the total area of Chile (in this raster, at least).

````@example cellarea
total_area = sum(skipmissing(masked_areas))
````

And we can convert that to an average by dividing by the total area:

````@example cellarea
avg_precip = total_precip / total_area
````

According to the internet, Chile gets about 100mm of rain per square meter in June, so our statistic seems pretty close.

Let's see what happens if we don't account for cell areas.  An equivalent assumption would be that all cells have the same area.

````@example cellarea
bad_total_precip = sum(skipmissing(masked_precip))
bad_avg_precip = bad_total_precip / length(collect(skipmissing(masked_precip)))
````

This is misestimated!  This is why it's important to account for cell areas when computing averages.

!!! note
    If you made it this far, congratulations!

    It's interesting to note that we've replicated the workflow of `zonal` here.
    `zonal` is a more general function that can be used to compute any function over geometries,
    and it has multithreading built in.

    But fundamentally, this is all that `zonal` is doing under the hood -
    masking and cropping the raster to the geometry, and then computing the statistic.

## Summary

In this tutorial, we've seen how to compute the spatial mean of a raster, and how to account for the fact that raster cells do not always have the same area.

We've also seen how to use the `cellarea` function to compute the area of each cell in a raster, and how to use the `mask` function to get the data within a geometry.

We've seen that the spatial mean is not the same as the arithmetic mean, and that we need to account for the area of each cell when computing the average.

