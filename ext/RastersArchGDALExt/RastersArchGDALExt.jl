module RastersArchGDALExt

import ArchGDAL,
    DiskArrays,
    Extents,
    Missings

using DimensionalData,
    GeoFormatTypes,
    GeoInterface,
    Rasters

using Rasters.LookupArrays
using Rasters.Dimensions
using Rasters: GDALsource, AbstractProjected, RasterStackOrArray, FileArray,
    RES_KEYWORD, SIZE_KEYWORD, CRS_KEYWORD, EXPERIMENTAL, GDAL_EMPTY_TRANSFORM, 
    GDAL_TOPLEFT_X, GDAL_WE_RES, GDAL_ROT1, GDAL_TOPLEFT_Y, GDAL_ROT2, GDAL_NS_RES

import Rasters: reproject, resample, warp

const RA = Rasters
const DD = DimensionalData
const DA = DiskArrays
const GI = GeoInterface
const LA = LookupArrays

include("gdal_source.jl")
include("reproject.jl")
include("resample.jl")
include("warp.jl")

end
