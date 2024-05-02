module RastersZarrDatasetsExt

using Rasters, CommonDataModel
using ZarrDatasets: ZarrDatasets as ZD

import DiskArrays,
    FillArrays,
    Extents,
    GeoInterface,
    Missings

using Dates, 
    DimensionalData,
    GeoFormatTypes

using Rasters.Lookups
using Rasters.Dimensions
using Rasters: Zarrsource

using ZarrDatasets
using CommonDataModel: AbstractDataset

const RA = Rasters
const DD = DimensionalData
const DA = DiskArrays
const GI = GeoInterface
const LA = Lookups

include("zarrdatasets_source.jl")

end
