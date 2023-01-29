# [Plots.jl](https://github.com/JuliaPlots/Plots.jl) is fully supported by
# Rasters.jl, with recipes for plotting `Raster` and `RasterStack` provided. `plot`
# will plot a heatmap with axes matching dimension values. If `mappedcrs` is used,
# converted values will be shown on axes instead of the underlying `crs` values.
# `contourf` will similarly plot a filled contour plot.

# Pixel resolution is limited to allow loading very large files quickly. `max_res` 
# specifies the maximum pixel resolution to show on the longest axis of the array.
# It can be set manually to change the resolution (e.g. for large or high-quality plots):
using Rasters
const RS = Rasters
using CairoMakie
CairoMakie.activate!()

A = Raster(WorldClim{BioClim}, 5)
rplot(A)

# ## Loading and plotting data

# Our first example simply loads a file from disk and plots it.

# This netcdf file only has one layer, if it has more we could use `RasterStack`
# instead.

url = "https://www.unidata.ucar.edu/software/netcdf/examples/tos_O1_2001-2002.nc";
filename = download(url, "tos_O1_2001-2002.nc");
A = Raster(filename)

# Objects with Dimensions other than `X` and `Y` will produce multi-pane plots.
# Here we plot every third month in the first year in one plot:

rasd = A[Ti=1:3:12]

rplot(rasd; colormap= :seaborn_icefire_gradient)

# Now plot the ocean temperatures around the Americas in the first month of 2001.
# Notice we are using lat/lon coordinates and date/time instead of regular
# indices. The time dimension uses `DateTime360Day`, so we need to load CFTime.jl
# to index it with `Near`.

using CFTime
tband = A[Ti(Near(DateTime360Day(2001, 01, 17))), Y(-60.0 .. 90.0), X(45.0 .. 190.0)]
rplot(tband)

# Other plot functions and sliced objects that have only one `X`/`Y`/`Z` dimension
# fall back to generic DimensionalData.jl plotting, which will still correctly
# label plot axes.

# ## Plot a contourf plot

using Statistics
## Take the mean
mean_tos = mean(A; dims=Ti)
mean_tos = replace_missing(mean_tos[Ti(1)], NaN)
x, y = lookup(mean_tos, X), lookup(mean_tos, Y)

contourf(x, y, mean_tos.data;
    axis = (; aspect = DataAspect()),
    )

# Write the mean values to disk

# ```
# write("mean_tos.nc", mean_tos)
# ```

# Plotting recipes in DimensionalData.jl are the fallback for Rasters.jl when the
# object doesn't have 2 `X`/`Y`/`Z` dimensions, or a non-spatial plot command is
# used. So (as a random example) we could plot a transect of ocean surface
# temperature at 20 degree latitude:

dnear = A[Y(Near(20.0)), Ti(1)]
x = lookup(dnear, X)
lines(x, dnear.data)