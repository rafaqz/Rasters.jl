module RastersZarrDatasetsExt

using Rasters
using ZarrDatasets

using ZarrDatasets: ZarrDatasets as ZD
using Rasters: Zarrsource

const RA = Rasters

Rasters.is_loaded(::Rasters.ZarrDatasetsExt) = true

include("zarrdatasets_source.jl")

end
