module RastersZarrDatasetsExt

using Rasters
using ZarrDatasets

using ZarrDatasets: ZarrDatasets as ZD
using Rasters: Zarrsource

const RA = Rasters

include("zarrdatasets_source.jl")

end
