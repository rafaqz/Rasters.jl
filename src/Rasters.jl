module Rasters

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end Rasters

using Dates

# Load first to fix StaticArrays invalidations
import CoordinateTransformations
import DimensionalData

import Adapt,
       ArchGDAL,
       ColorTypes,
       ConstructionBase,
       DiskArrays,
       Extents,
       FillArrays,
       Flatten,
       GeoInterface,
       HDF5,
       OffsetArrays,
       ProgressMeter,
       MakieCore,
       Missings,
       Mmap,
       NCDatasets,
       RecipesBase,
       Reexport,
       Setfield

# This symbol is only defined on Julia versions that support extensions.
@static if !isdefined(Base, :get_extension)
    using Requires
end

Reexport.@reexport using DimensionalData, GeoFormatTypes, RasterDataSources

using DimensionalData.Tables,
      DimensionalData.LookupArrays,
      DimensionalData.Dimensions
      DimensionalData.LookupArrays.IntervalSets

using DimensionalData: Name, NoName
using .Dimensions: StandardIndices, DimTuple
using .LookupArrays: LookupArrayTuple 

using RecipesBase: @recipe, @series
using Base: tail, @propagate_inbounds

using Setfield: @set, @set!
using ColorTypes: RGB

export AbstractRaster, Raster
export AbstractRasterStack, RasterStack
export AbstractRasterSeries, RasterSeries
export Projected, Mapped
export Band
export missingval, boolmask, missingmask, replace_missing, replace_missing!,
       aggregate, aggregate!, disaggregate, disaggregate!, mask, mask!, 
       resample, warp, zonal, crop, extend, trim, slice, combine, points,
       classify, classify!, mosaic, mosaic!, extract, rasterize, rasterize!,
       coverage, coverage!, setcrs, setmappedcrs
export crs, mappedcrs, mappedindex, mappedbounds, projectedindex, projectedbounds
export reproject, convertlookup


const DD = DimensionalData
const DA = DiskArrays
const GI = GeoInterface
const LA = LookupArrays

# DimensionalData documentation urls
const DDdocs = "https://rafaqz.github.io/DimensionalData.jl/stable/api"
const DDdimdocs = joinpath(DDdocs, "#DimensionalData.Dimension")
const DDarraydocs = joinpath(DDdocs, "#DimensionalData.AbstractDimensionalArray")
const DDabssampleddocs = joinpath(DDdocs, "#DimensionalData.AbstractSampled")
const DDsampleddocs = joinpath(DDdocs, "#DimensionalData.Sampled")
const DDlocusdocs = joinpath(DDdocs, "#DimensionalData.Locus")
const DDselectordocs = joinpath(DDdocs, "#DimensionalData.Selector")
const DDtidocs = joinpath(DDdocs, "#DimensionalData.Ti")

const EXPERIMENTAL = """
    WARNING: This feature is experimental. It may change in future versions, and may
    not be 100% reliable in all cases. Please file github issues if problems occur.
    """

include("lookup.jl")
include("dimensions.jl")
include("filearray.jl")
include("filestack.jl")
include("openstack.jl")
include("array.jl")
include("stack.jl")
include("series.jl")

const RasterStackOrArray = Union{AbstractRasterStack,AbstractRaster}
const RasterSeriesOrStack = Union{AbstractRasterSeries,AbstractRasterStack}

include("utils.jl")
include("skipmissing.jl")
include("polygon_ops.jl")
include("table_ops.jl")
include("create.jl")
include("read.jl")
include("write.jl")
include("show.jl")
include("plotrecipes.jl")
include("sectorlock.jl")


include("methods/mask.jl")
include("methods/rasterize.jl")
include("methods/aggregate.jl")
include("methods/classify.jl")
include("methods/crop_extend.jl")
include("methods/coverage.jl")
include("methods/extract.jl")
include("methods/mosaic.jl")
include("methods/points.jl")
include("methods/replace_missing.jl")
include("methods/reproject.jl")
include("methods/resample.jl")
include("methods/slice_combine.jl")
include("methods/trim.jl")
include("methods/warp.jl")
include("methods/zonal.jl")

include("sources/sources.jl")
include("sources/grd.jl")
include("sources/smap.jl")
include("sources/ncdatasets.jl")
include("sources/gdal.jl")
include("sources/rasterdatasources.jl")

# extensions

function __init__()
    @static if !isdefined(Base, :get_extension)
        @require ArchGDAL = "" include("../ext/RastersArchGDALExt.jl")
        @require CoordinateTransformations = "" include("../ext/RastersCoordinateTransformationsExt.jl")
        @require Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("../ext/RastersMakie.jl")
        @require NCDatasets = "" include("../ext/RastersNCDatasetsExt.jl")
        @require RasterDataSources = "" include("../ext/RastersRasterDataSources.jl")
    end
end

end
