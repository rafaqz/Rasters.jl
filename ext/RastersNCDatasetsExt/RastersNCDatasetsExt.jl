module RastersNCDatasetsExt

@static if isdefined(Base, :get_extension) # julia < 1.9
    using Rasters, NCDatasets
else    
    using ..Rasters, ..NCDatasets
end

import DiskArrays,
    FillArrays,
    Extents,
    GeoInterface,
    Missings,

using Dates, 
    DimensionalData,
    GeoFormatTypes

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
