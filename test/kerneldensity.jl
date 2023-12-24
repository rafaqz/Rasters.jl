using Test, Rasters, KernelDensity, Extents

points = collect(zip(rand(100) .* 10, rand(100)))
# Using just res, detecting extent from points
res_rast = kerneldensity(points; res=(0.1, 0.01))
# This doesn't exactly match, someting is wrong with `_extent2dims`
@test_broken map(step, dims(res_rast)) == (0.1, 0.01)
# Using just size, detecting extent from points
size_rast = kerneldensity(points; size=250)
@test size(size_rast) == (250, 250)
# Using a defined extent with res
to = Extents.Extent(X=(0, 10), Y=(0, 5))
ext_res_rast = kerneldensity(points; to, res=0.05)
@test size(ext_rast) == (200, 100)
@test map(step, dims(ext_rast)) == (0.05, 0.05)
# Using a defined extent with size
to = Extents.Extent(X=(0, 10), Y=(0, 5))
ext_res_rast = kerneldensity(points; to, size=(200, 100))
@test size(ext_rast) == (200, 100)
@test map(step, dims(ext_rast)) == (0.05, 0.05)
# Using a defined DimArray
to = rand(X(0:0.01:10), Y(0:0.001:1))
template_rast = kerneldensity(points; to)
@test size(to) == size(template_rast)
@test map(step, dims(template_rast)) == (0.01, 0.001)

# Test plots
# using GLMakie
# p = Makie.heatmap(res_rast)
# Makie.scatter!(points)
# p = Makie.heatmap(size_rast)
# Makie.scatter!(points)
# p = Makie.heatmap(ext_rast)
# Makie.scatter!(points)
# p = Makie.heatmap(template_rast)
# Makie.scatter!(points)
