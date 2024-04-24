module RastersZarrDatasetsExt

@static if isdefined(Base, :get_extension) # julia < 1.9
    using Rasters, ZarrDatasets, CommonDataModel
else    
    using ..Rasters, ..GRIBDatasets, ..CommonDataModel
end

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
