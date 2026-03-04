module RastersRasterDataSourcesExt

using Rasters, RasterDataSources

using Rasters.Lookups
using Rasters.Dimensions
using Rasters: DiskArrays, FillArrays

const RA = Rasters
const RDS = RasterDataSources
const DA = DiskArrays

Rasters.is_loaded(::Rasters.RasterDataSourcesExt) = true

include("constructors.jl")

end # Module
