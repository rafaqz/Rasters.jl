# ## Plotting in Makie

# Plotting in Makie works somewhat differently than Plots, since the recipe system is different.
# You can pass a 2-D raster to any surface-like function (heatmap, contour, contourf,
# or even surface for a 3D plot) with ease.

# ### 2-D rasters in Makie

using CairoMakie, Makie
using Rasters, RasterDataSources, ArchGDAL
A = Raster(WorldClim{BioClim}, 5) # this is a 3D raster, so is not accepted.

fig = Figure()
plot(fig[1, 1], A)
contour(fig[1, 2], A)
ax = Axis(fig[2, 1]; aspect = DataAspect())
contourf!(ax, A)
surface(fig[2, 2], A) # even a 3D plot works!
fig

# ### 3-D rasters in Makie

# !!! warning
#       This interface is experimental, and unexported for that reason.  It may break at any time!

# Just as in Plots, 3D rasters are treated as a series of 2D rasters, which are tiled and plotted.  

# You can use `Rasters.rplot` to visualize 3D rasters or RasterStacks in this way.  An example is below:

stack = RasterStack(WorldClim{Climate}; month = 1)
Rasters.rplot(stack; Axis = (aspect = DataAspect(),),)


# You can pass any theming keywords in, which are interpreted by Makie appropriately.

# The plots seem a little squished here.  We provide a Makie theme which makes text a little smaller
# and has some other space-efficient attributes:

Makie.set_theme!(Rasters.theme_rasters())
Rasters.rplot(stack)


# reset theme
Makie.set_theme!() 

# ### Plotting with `Observable`s, animations

# `Rasters.rplot` should support Observable input out of the box, but the dimensions of that input
# must remain the same - i.e., the element names of a RasterStack must remain the same.

Makie.set_theme!(Rasters.theme_rasters())

# `stack` is the WorldClim climate data for January
stack_obs = Observable(stack)
fig = Rasters.rplot(stack_obs;
    Colorbar=(; height=Relative(0.75), width=5)
) 
record(fig, "rplot.mp4", 1:12; framerate = 3) do i
    stack_obs[] = RasterStack(WorldClim{Climate}; month = i)
end 

# ![](rplot.mp4)

Makie.set_theme!() # reset theme

# ```@docs
# Rasters.rplot
# ```

# ## Using vanilla Makie

using Rasters, RasterDataSources

# The data
layers = (:evenness, :range, :contrast, :correlation)
st = RasterStack(EarthEnv{HabitatHeterogeneity}, layers)
ausbounds = X(100 .. 160), Y(-50 .. -10) # Roughly cut out australia
aus = st[ausbounds...] |> Rasters.trim

# The plot
# colorbar attributes
colormap = :batlow
flipaxis = false
tickalign=1
width = 13
ticksize = 13
# figure
with_theme(theme_ggplot2()) do 
    fig = Figure(resolution=(800, 600))
    axs = [Axis(fig[i,j], xlabel = "lon", ylabel = "lat") for i in 1:2 for j in 1:2]
    for (i, l) in enumerate(layers) 
        Makie.heatmap!(axs[i], aus[l]; colormap)
        axs[i].title = string(l) 
    end
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

# ![aus_trim](aus_trim.png)
