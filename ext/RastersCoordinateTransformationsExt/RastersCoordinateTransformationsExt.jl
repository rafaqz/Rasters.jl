module RastersCoordinateTransformationsExt

@static if isdefined(Base, :get_extension) # julia < 1.9
    using Rasters, CoordinateTransformations
else    
    using ..Rasters, ..CoordinateTransformations
end

using DimensionalData
using Rasters.LookupArrays
using Rasters.Dimensions

import Rasters: AffineProjected, GDAL_EMPTY_TRANSFORM, GDAL_TOPLEFT_X, 
                GDAL_WE_RES, GDAL_ROT1, GDAL_TOPLEFT_Y, GDAL_ROT2, GDAL_NS_RES
const RA = Rasters
const DD = DimensionalData
const LA = LookupArrays


include("affineprojected.jl")

end # module
