module RastersRasterDataSourcesExt

using Rasters, RasterDataSources

using Rasters.Lookups
using Rasters.Dimensions

const RA = Rasters
const RDS = RasterDataSources

Rasters.is_loaded(::Rasters.RasterDataSourcesExt) = true

include("constructors.jl")

end # Module
