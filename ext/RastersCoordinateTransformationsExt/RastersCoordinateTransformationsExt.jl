module RastersCoordinateTransformationsExt

using DimensionalData,
    CoordinateTransformations,
    Rasters

using Rasters.LookupArrays
using Rasters.Dimensions

const RA = Rasters
const DD = DimensionalData
const LA = LookupArrays


include("affineprojected.jl")
include("geotransform.jl")

end # module
