module RastersGRIBDatasetsExt

@static if isdefined(Base, :get_extension) # julia < 1.9
    using Rasters, GRIBDatasets, CommonDataModel
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

using Rasters.LookupArrays
using Rasters.Dimensions
using Rasters: GRIBsource

using CommonDataModel: AbstractDataset

const RA = Rasters
const DD = DimensionalData
const DA = DiskArrays
const GI = GeoInterface
const LA = LookupArrays

include("gribdatasets_source.jl")

end
