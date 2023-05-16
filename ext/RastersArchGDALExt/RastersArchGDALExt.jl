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
    RES_KEYWORD, SIZE_KEYWORD, CRS_KEYWORD, EXPERIMENTAL

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
