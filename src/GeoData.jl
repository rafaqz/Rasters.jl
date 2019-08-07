module GeoData

using Reexport, RecipesBase
@reexport using CoordinateReferenceSystemsBase, DimensionalData

using DimensionalData: formatdims, slicedims, val
import DimensionalData: dims, refdims, rebuild

export AbstractGeoArray, GeoArray
export AbstractGeoData, GeoData
export AbstractGeoStack, GeoStack
export coordinates, coordinates!, missingval, metadata, replace_missing


include("types.jl")
include("interface.jl")
include("coordinates.jl")
include("geoarray.jl")
include("plotrecipes.jl")

end
