module RastersNCDatasetsExt

import DiskArrays,
    Extents,
    Missings,
    NCDatasets

using Dates, 
    DimensionalData,
    GeoFormatTypes,
    GeoInterface,
    Rasters

using Rasters.LookupArrays
using Rasters.Dimensions
using Rasters: NCDsource, AbstractProjected, RasterStackOrArray, FileArray, FileStack, OpenStack, DimTuple, Key,
    RES_KEYWORD, SIZE_KEYWORD, CRS_KEYWORD, EXPERIMENTAL

import Rasters: Raster, create, crs, dims, refdims, metadata, missingval, reproject, resample, warp,
    _open, _writeable_missing, _metadatadict, cleanreturn

const RA = Rasters
const DD = DimensionalData
const DA = DiskArrays
const GI = GeoInterface
const LA = LookupArrays

include("ncdatasets_source.jl")

end
