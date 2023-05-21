module RastersRasterDataSourcesExt

using Rasters, RasterDataSources

using RasterDataSources: RasterDataSource
using Rasters.LookupArrays
using Rasters.Dimensions

const RA = Rasters
const RDS = RasterDataSources

include("constructors.jl")

end # Module
