module RastersGRIBDatasetsExt

using Rasters
using GRIBDatasets

using Rasters: GRIBsource

const RA = Rasters
const GDS = GRIBDatasets

include("gribdatasets_source.jl")

end
