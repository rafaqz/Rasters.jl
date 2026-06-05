# Rasters with Makie.jl

## Setup

Install the required packages by entering the Julia REPL package mode (press `]`) and typing:

```julia
add Rasters CairoMakie Makie RasterDataSources ArchGDAL
```

or from a script/notebook: 

````julia
using Pkg
Pkg.add(["Rasters", "CairoMakie", "Makie", "RasterDataSources", "ArchGDAL"])
````

To download data you will need to specify a folder to put it in. You can do this by assigning the environment variable RASTERDATASOURCES_PATH: 

````julia
ENV["RASTERDATASOURCES_PATH"] = "/home/user/Data/" # your path here
````

## Plotting in Makie

Plotting in Makie works somewhat differently than Plots, since the recipe system is different.
You can pass a 2-D raster to any surface-like function (heatmap, contour, contourf,
or even surface for a 3D plot) with ease.

## 2-D rasters in Makie

We'll start with a single 2-D raster: BioClim variable 5 from the WorldClim dataset, which is the maximum temperature of the warmest month (in °C) on a global grid. 


````@example makie
using CairoMakie, Makie
using Rasters, RasterDataSources, ArchGDAL
A = Raster(WorldClim{BioClim}, 5)
````

````@example makie
fig, ax, _ = plot(A)
contour(fig[1, 2], A)
ax = Axis(fig[2, 1]; aspect = DataAspect())
contourf!(ax, A)
surface(fig[2, 2], A; axis = (type = LScene,)) # even a 3D plot works!
fig
````

## 3-D rasters in Makie

!!! warning
      This interface is experimental, and unexported for that reason.  It may break at any time!

Just as in Plots, 3D rasters are treated as a series of 2D rasters, which are tiled and plotted.  

You can use `Rasters.rplot` to visualize 3D rasters or RasterStacks in this way.  An example is below:

````@example makie
stack = RasterStack(WorldClim{Climate}; month = 1)
Rasters.rplot(stack; Axis = (aspect = DataAspect(),),)
````

You can pass any theming keywords in, which are interpreted by Makie appropriately.

### Plotting with `Observable`s, animations

An `Observable` is a container object whose stored value you can update interactively. You can create functions or other observables that are executed whenever an observable changes. This is how Makie supports interactive and animated plots. 

`Rasters.rplot` should support Observable input out of the box, but the dimensions of that input must remain the same - i.e., the element names of a RasterStack must remain the same.

````@example makie
# `stack` is the WorldClim climate data for January
stack_obs = Observable(stack)
fig = Rasters.rplot(stack_obs;
    Colorbar=(; height=Relative(0.75), width=5)
) 
record(fig, "rplot.mp4", 1:12; framerate = 3) do i
    stack_obs[] = RasterStack(WorldClim{Climate}; month = i)
end
````

```@raw html
<!-- <video src="./rplot.mp4" controls="controls" autoplay="autoplay"></video> -->
```

````@example makie
Makie.set_theme!() # reset theme
````

```@docs
Rasters.rplot
```

## Using vanilla Makie

So far we've leaned on Rasters' Makie recipes, which take a raster directly and set up the axes and layout for us. For full control over the figure, we can instead build the scaffolding manually.

````@example makie
using Rasters, RasterDataSources
````
The data:

````@example makie
layers = (:evenness, :range, :contrast, :correlation) # tuple of the four variable names from HabitatHeterogeneity we want to use
st = RasterStack(EarthEnv{HabitatHeterogeneity}, layers) # load the four layers together as a RasterStack
ausbounds = X(100 .. 160), Y(-50 .. -10) # Roughly cut out Australia using the .. selector from DimensionalData.jl
aus = st[ausbounds...] |> Rasters.trim # crop to the bounds, then trim away the empty (missing) edge rows/columns
````
!!! note
      Rasters extends [`DimensionalData`](https://rafaqz.github.io/DimensionalData.jl/stable/)'s `DimArray` and `DimStack` to build `Raster` and `RasterStack` objects, and uses its selectors such as `..`

The plot:

````@example makie
# colorbar attributes
colormap = :batlow
flipaxis = false
tickalign=1
width = 13
ticksize = 13
# figure
with_theme(theme_dark()) do 
    fig = Figure(; size=(600, 600), backgroundcolor=:transparent)
    axs = [Axis(fig[i,j], xlabel = "lon", ylabel = "lat",
        backgroundcolor=:transparent) for i in 1:2 for j in 1:2]
    plt = [Makie.heatmap!(axs[i], aus[l]; colormap) for (i, l) in enumerate(layers)]
    for (i, l) in enumerate(layers) axs[i].title = string(l) end
    hidexdecorations!.(axs[1:2]; grid=false, ticks=false)
    hideydecorations!.(axs[[2,4]]; grid=false, ticks=false)
    Colorbar(fig[1, 0], plt[1]; flipaxis, tickalign, width, ticksize)
    Colorbar(fig[1, 3], plt[2]; tickalign, width, ticksize)
    Colorbar(fig[2, 0], plt[3]; flipaxis, tickalign, width, ticksize)
    Colorbar(fig[2, 3], plt[4]; tickalign, width, ticksize)
    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 5)
    Label(fig[0, :], "RasterStack of EarthEnv HabitatHeterogeneity layers, trimmed to Australia")
    fig 
end
save("aus_trim.png", current_figure());
````

![aus_trim](aus_trim.png)