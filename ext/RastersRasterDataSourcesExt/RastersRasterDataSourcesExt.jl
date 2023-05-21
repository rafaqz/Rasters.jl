module RastersRasterDataSourcesExt

@static if isdefined(Base, :get_extension) # julia < 1.9
    using Rasters, RasterDataSources
else    
    using ..Rasters, ..RasterDataSources
end

# using RasterDataSources: RasterDataSource
using Rasters.LookupArrays
using Rasters.Dimensions

const RA = Rasters
const RDS = RasterDataSources

include("constructors.jl")

end # Module
