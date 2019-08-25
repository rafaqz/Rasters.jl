module GeoData

using Mixers, RecipesBase, Reexport, Requires
@reexport using CoordinateReferenceSystemsBase, DimensionalData


using DimensionalData: Time, formatdims, slicedims, basetype, dims2indices

import DimensionalData: val, dims, refdims, metadata, rebuild, select, selectview
import CoordinateReferenceSystemsBase: crs

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
include("utils.jl")

DimensionalData.@dim Band

# function __init__()
    # @require HDF5="f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f" begin
        # This section is for sources that rely on HDF5
        # Not simply any HDF5.
        include("sources/smap.jl")
    # end
    # @require NCDatasets="85f8d34a-cbdd-5861-8df4-14fed0d494ab" 
     include("sources/ncdatasets.jl") 
    # @require ArchGDAL="c9ce4bd3-c3d5-55b8-8973-c0e20141b8c3" 
    include("sources/gdal.jl")
# end

end
