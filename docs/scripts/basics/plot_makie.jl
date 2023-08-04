# ## Plotting in Makie

# Plotting in Makie works somewhat differently than Plots, since the recipe system is different.
# You can pass a 2-D raster to any surface-like function (heatmap, contour, contourf,
# or even surface for a 3D plot) with ease.

# ### 2-D rasters in Makie

using CairoMakie, Makie
using Rasters, RasterDataSources, ArchGDAL
A = Raster(WorldClim{BioClim}, 5) # this is a 3D raster, so is not accepted.
B = A[:, :, 1] # this converts to a 2D raster which Makie accepts!

fig = Figure()
plot(fig[1, 1], B)
contour(fig[1, 2], B)
ax = Axis(fig[2, 1]; aspect = DataAspect())
contourf!(ax, B)
surface(fig[2, 2], B) # even a 3D plot works!
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

CairoMakie.set_theme!(Rasters.theme_rasters())
Rasters.rplot(stack)


# reset theme
using Makie
Makie.set_theme!() 

# ### Plotting with `Observable`s, animations

# `Rasters.rplot` should support Observable input out of the box, but the dimensions of that input
# must remain the same - i.e., the element names of a RasterStack must remain the same.

CairoMakie.set_theme!(Rasters.theme_rasters())

stack_obs = Observable(stack)
fig = Rasters.rplot(stack_obs;
    Colorbar = (; height= Relative(0.75), width=5,)) # `stack` is the WorldClim climate data for January
record(fig, "rplot.mp4", 1:12; framerate = 3) do i
    stack_obs[] = RasterStack(WorldClim{Climate}; month = i)
end 

# ![](rplot.mp4)

using Makie
Makie.set_theme!() # reset theme

# ```@docs
# Rasters.rplot
# ```