module RastersNCDatasetsExt

import DiskArrays,
    FillArrays,
    Extents,
    GeoInterface,
    Missings,
    NCDatasets

using Dates, 
    DimensionalData,
    GeoFormatTypes,
    Rasters

using Rasters.LookupArrays
using Rasters.Dimensions
using Rasters: NCDsource

const RA = Rasters
const DD = DimensionalData
const DA = DiskArrays
const GI = GeoInterface
const LA = LookupArrays

include("ncdatasets_source.jl")

end
