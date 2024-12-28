## `cellarea`

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

````@example cella
using Rasters, Proj
````

To construct a raster, we'll need to specify the x and y dimensions.  These are called "lookups" in Rasters.jl.
using Rasters.Lookups
We can now construct the x and y lookups.  Here we'll use a start-at-one, step-by-five grid.
Note that we're specifying that the "sampling", i.e., what the coordinates actually mean, 
is `Intervals(Start())`, meaning that the coordinates are the starting point of each interval.

This is in contrast to `Points()` sampling, where each index in the raster represents the value at a sampling point;
here, each index represents a grid cell, which is defined by the coordinate being at the start.

````@example cella
x = X(1:5:30; sampling = Intervals(Start()), crs = EPSG(4326))
y = Y(50:5:80; sampling = Intervals(Start()), crs = EPSG(4326))
````

I have chosen the y-range here specifically so we can show the difference between spherical and planar `cellarea`.

````@example cella
ras = Raster(ones(x, y); crs = EPSG(4326))
````

We can just call `cellarea` on this raster, which returns cell areas in meters, on Earth, assuming it's a sphere:

````@example cella
cellarea(ras)
````

and if we plot it, you can see the difference in cell area as we go from the equator to the poles:

````@example cella
using Makie, CairoMakie# , GeoMakie
heatmap(ras; axis = (; aspect = DataAspect()))
````

We can also try this using the planar method, which simply computes the area of the rectangle using `area = x_side_length * y_side_length`:

````@example cella
cellarea(ras, Planar())
````

Note that this is of course wildly inaccurate for a geographic dataset - but if you're working in a projected coordinate system, like polar stereographic or Mercator, this can be very useful (and a _lot_ faster)!