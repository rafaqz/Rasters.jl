# Rasterize with crazy things

A benefit of being written in pure Julia is that you can do anything you want with 
`rasterize` - you can rasterize to any type or combination of types that you can get to work.

Let's take an example.  Say you want to rasterize all the countries in the world, but 
where two countries overlap, you want to save the indices of both countries.

A simple but inefficient way to do this is to have each pixel of the raster actually be an
array of integers, which stores the index of each geometry that touches the pixel.  Then you
can control exactly what you want to do.

First, let's get some data.  This is a feature collection of all countries in the world.

```@example crazy
using NaturalEarth: naturalearth
using CairoMakie # plotting

countries = naturalearth("admin_0_countries", 10)
poly(countries.geometry; strokewidth = 1)
```

Next, we can rasterize it.  Note that if you're going beyond the standard datatypes like `Int`,
`Float64`, etc., you will need to specify some extra kwargs to `rasterize`, so Rasters doesn't
have to guess what you want. 

```@example crazy
using Rasters

ras = rasterize(
    countries;
    op = vcat, 
    fill = [[i] for i in 1:length(countries.geometry)],
    boundary = :touches, # this one is just for the plot to look good
    res = 0.5, # half-degree resolutoin
    # Below are the kwargs you have to provide if using a
    # custom type that might not have e.g. `zero` defined
    eltype = Vector{Int},
    missingval = Int[], # empty vector of Int
    init = Int[],       # empty vector of Int
    progress = false, # hide
)
```

You can also see that this has a _lot_ of unique values:
```@example crazy
length(unique(ras))
```

It's clearest if you plot a heatmap of the lengths of each cell:

```@example crazy
heatmap(length.(ras))
```