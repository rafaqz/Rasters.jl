module RastersProjExt

using Rasters, Proj

import DiskArrays,
    Extents,
    Missings

using DimensionalData,
    GeoFormatTypes,
    GeoInterface


using Rasters.Lookups
using Rasters.Dimensions
using Rasters: GDALsource, AbstractProjected, RasterStackOrArray, FileArray, NoKW,
    RES_KEYWORD, SIZE_KEYWORD, CRS_KEYWORD, FILENAME_KEYWORD, SUFFIX_KEYWORD, EXPERIMENTAL,
    GDAL_EMPTY_TRANSFORM, GDAL_TOPLEFT_X, GDAL_WE_RES, GDAL_ROT1, GDAL_TOPLEFT_Y, GDAL_ROT2, GDAL_NS_RES,
    _no_crs_error

import Rasters: reproject, _spherical_cellarea, nokw

import LinearAlgebra: dot, cross

const RA = Rasters
const DD = DimensionalData
const DA = DiskArrays
const GI = GeoInterface
const LA = Lookups

include("cellarea.jl")
include("reproject.jl")
end