module RastersRDSZarrExt

using Rasters
using RasterDataSources
using ZarrDatasets
using Zarr

using Rasters: Zarrsource
using Rasters.Lookups
using Rasters.Dimensions
using ZarrDatasets: ZarrDatasets as ZD

const RA = Rasters
const RDS = RasterDataSources

include("cachedcloudsource.jl")

end
