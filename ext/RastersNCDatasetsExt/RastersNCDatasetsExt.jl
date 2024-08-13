module RastersNCDatasetsExt

using Rasters
using NCDatasets
using CommonDataModel
using Dates 
using DimensionalData

import Missings

using Rasters.Lookups
using Rasters.Dimensions
using Rasters: CDMsource, NCDsource, NoKW, nokw, isnokw

using CommonDataModel: AbstractDataset

const NCD = NCDatasets
const CDM = CommonDataModel
const RA = Rasters
const DD = DimensionalData
const LA = Lookups

include("ncdatasets_source.jl")

end
