module RastersMakieExt

using Makie
using Rasters

using Rasters.DimensionalData
using Rasters.Dimensions

Rasters.is_loaded(::Rasters.MakieExt) = true

include("plotrecipes.jl")

end
