module RastersGRIBDatasetsExt

using Rasters
using GRIBDatasets
using DiskArrays

using Rasters: GRIBsource

const RA = Rasters
const DA = DiskArrays
const GDS = GRIBDatasets

Rasters.is_loaded(::Type{Rasters.GRIBDatasetsExt}) = true

include("gribdatasets_source.jl")

end
