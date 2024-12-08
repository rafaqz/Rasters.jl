module RastersCoordinateTransformationsExt

using Rasters
using CoordinateTransformations
using DimensionalData
using Rasters.Lookups
using Rasters.Dimensions

import Rasters: AffineProjected, GDAL_EMPTY_TRANSFORM, GDAL_TOPLEFT_X, 
                GDAL_WE_RES, GDAL_ROT1, GDAL_TOPLEFT_Y, GDAL_ROT2, GDAL_NS_RES
const RA = Rasters
const DD = DimensionalData
const LA = Lookups

include("affineprojected.jl")

end # module
