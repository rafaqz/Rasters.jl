module GeoData

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end GeoData

using Missings,
      Mixers,
      RecipesBase,
      Reexport,
      Requires,
      Dates,
      Mmap

@reexport using DimensionalData, GeoFormatTypes

const DD = DimensionalData

using Base: tail
using DimensionalData: Forward, Reverse, formatdims, slicedims, basetypeof,
      dims2indices, indexorder, arrayorder, relationorder, StandardIndices, isrev

import DimensionalData: val, data, dims, refdims, metadata, rebuild, rebuildsliced,
                        name, label, units, bounds, sel2indices, mode, 
                        order, locus, span, sampling, forwardorder

export Metadata, ArrayMetadata, DimMetadata

export AbstractGeoArray, GeoArray

export AbstractGeoStack, GeoStack

export AbstractGeoSeries, GeoSeries

export Projected

export missingval, boolmask, missingmask, replace_missing, aggregate, aggregate!, crs, usercrs

export Lon, Lat, Vert, Band

@dim Lon XDim "Longitude" "Lon"
@dim Lat YDim "Latitude" "Lat"
@dim Vert ZDim "Vertical" "Vert"
@dim Band

include("interface.jl")
include("metadata.jl")
include("array.jl")
include("stack.jl")
include("series.jl")
include("utils.jl")
include("aggregate.jl")
include("methods.jl")
include("mode.jl")
include("sources/grd.jl")
include("plotrecipes.jl")


function __init__()
    @require HDF5="f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f" begin
        # This section is for sources that rely on HDF5, not simply any HDF5.
        include("sources/smap.jl")
    end
    @require NCDatasets="85f8d34a-cbdd-5861-8df4-14fed0d494ab" begin
        include("sources/ncdatasets.jl")
    end
    @require ArchGDAL="c9ce4bd3-c3d5-55b8-8973-c0e20141b8c3" begin
        include("reproject.jl")
        include("sources/gdal.jl")
    end
end

end
