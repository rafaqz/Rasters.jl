module GeoData

# Use the README as the module docs
@doc let 
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end GeoData

using Mixers, RecipesBase, Reexport, Requires, GeoFormatTypes, Dates

@reexport using DimensionalData


using DimensionalData: Time, X, Y, Z, Forward, Reverse, formatdims, slicedims, 
      basetypeof, dims2indices, @dim, indexorder, arrayorder
using Base: tail

import DimensionalData: val, dims, refdims, metadata, rebuild, rebuildsliced, name, label, units

export AbstractGeoArray, GeoArray
export AbstractGeoStack, GeoStack
export AbstractGeoSeries, GeoSeries
export missingval, mask, replace_missing
export Lon, Lat, Vert, Band

@dim Lon "Longitude"
@dim Lat "Latitude"
@dim Vert "Vertical"
@dim Band

include("interface.jl")
include("array.jl")
include("stack.jl")
include("series.jl")
include("plotrecipes.jl")
include("utils.jl")

function __init__()
    @require HDF5="f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f" begin
        # This section is for sources that rely on HDF5
        # Not simply any HDF5.
        include("sources/smap.jl")
    end
    @require NCDatasets="85f8d34a-cbdd-5861-8df4-14fed0d494ab" begin
        include("sources/ncdatasets.jl") 
    end
    @require ArchGDAL="c9ce4bd3-c3d5-55b8-8973-c0e20141b8c3" begin
        include("sources/gdal.jl")
    end
end

end
