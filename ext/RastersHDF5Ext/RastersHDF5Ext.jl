module RastersHDF5Ext

import DiskArrays,
    Extents,
    HDF5,
    Missings

using Dates, 
    DimensionalData,
    GeoFormatTypes,
    GeoInterface,
    Rasters

using Rasters.LookupArrays
using Rasters.Dimensions
using Rasters: SMAPsource, AbstractProjected, RasterStackOrArray, 
    FileArray, FileStack, OpenStack, DimTuple, Key, SMAPsource, RasterDiskArray,
    RES_KEYWORD, SIZE_KEYWORD, CRS_KEYWORD, EXPERIMENTAL, cleankeys

export smapseries

const RA = Rasters
const DD = DimensionalData
const DA = DiskArrays
const GI = GeoInterface
const LA = LookupArrays

include("smap_source.jl")

end
