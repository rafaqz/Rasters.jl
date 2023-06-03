module RastersHDF5Ext

@static if isdefined(Base, :get_extension) # julia < 1.9
    using Rasters, HDF5
else    
    using ..Rasters, ..HDF5
end

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
using Rasters: SMAPsource

export smapseries

const RA = Rasters
const DD = DimensionalData
const DA = DiskArrays
const GI = GeoInterface
const LA = LookupArrays

include("smap_source.jl")

end
