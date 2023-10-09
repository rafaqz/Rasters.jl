module RastersMakieExt

@static if isdefined(Base, :get_extension) # julia < 1.9
    using Makie, Rasters
else    
    using ..Makie, ..Rasters
end

using Rasters.DimensionalData
using Rasters.Dimensions

include("plotrecipes.jl")

end
