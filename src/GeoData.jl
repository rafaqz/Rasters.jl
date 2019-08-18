module GeoData

using Reexport, RecipesBase
@reexport using CoordinateReferenceSystemsBase, DimensionalData

using DimensionalData: formatdims, slicedims, basetype, val
import DimensionalData: dims, refdims, rebuild

export AbstractGeoArray, GeoArray
export AbstractGeoStack, GeoStack
export AbstractGeoSeries, GeoSeries
export coordinates, coordinates!, missingval, metadata, replace_missing


include("interface.jl")
include("array.jl")
include("stack.jl")
include("series.jl")
include("coordinates.jl")
include("plotrecipes.jl")

end
