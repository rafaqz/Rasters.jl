module RastersArchGDALExt

@static if isdefined(Base, :get_extension) # julia < 1.9
    using Rasters, ArchGDAL, CommonDataModel
else    
    using ..Rasters, ..ArchGDAL, ..CommonDataModel
end

import DiskArrays,
    Extents,
    Missings

using DimensionalData,
    GeoFormatTypes,
    GeoInterface

using Rasters.Lookups
using Rasters.Dimensions
using Rasters: GDALsource, AbstractProjected, AbstractRaster, AbstractRasterStack, 
    RasterStackOrArray, FileArray, NoKW,
    RES_KEYWORD, SIZE_KEYWORD, CRS_KEYWORD, FILENAME_KEYWORD, SUFFIX_KEYWORD, EXPERIMENTAL,
    GDAL_EMPTY_TRANSFORM, GDAL_TOPLEFT_X, GDAL_WE_RES, GDAL_ROT1, GDAL_TOPLEFT_Y, GDAL_ROT2, GDAL_NS_RES,
    _no_crs_error

import Rasters: reproject, resample, warp, cellsize, nokw

const RA = Rasters
const DD = DimensionalData
const DA = DiskArrays
const GI = GeoInterface
const LA = Lookups
const CDM = CommonDataModel

include("cellsize.jl")
include("gdal_source.jl")
include("reproject.jl")
include("resample.jl")
include("warp.jl")

end
