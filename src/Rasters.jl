module Rasters

# Use the README as the module docs
@doc let
    # path = joinpath(dirname(@__DIR__), "README.md")
    # include_dependency(path)
    # read(path, String)
end Rasters

using Dates

# Load first to fix StaticArrays invalidations
import DimensionalData

import Adapt,
       ColorTypes,
       CommonDataModel,
       ConstructionBase,
       DiskArrays,
       Extents,
       FillArrays,
       Flatten,
       GeoInterface,
       OffsetArrays,
       ProgressMeter,
       Missings,
       Mmap,
       RecipesBase,
       Reexport,
       Setfield

Reexport.@reexport using DimensionalData, GeoFormatTypes

using DimensionalData.Tables,
      DimensionalData.Lookups,
      DimensionalData.Dimensions
      DimensionalData.Lookups.IntervalSets

using DimensionalData: Name, NoName
using .Dimensions: StandardIndices, DimTuple
using .Lookups: LookupTuple

using RecipesBase: @recipe, @series
using Base: tail, @propagate_inbounds

import GeoInterface: crs

using Setfield: @set, @set!
using ColorTypes: RGB

using CommonDataModel: AbstractDataset, AbstractVariable

using DiskArrays: @implement_diskarray

export AbstractRaster, Raster
export AbstractRasterStack, RasterStack
export AbstractRasterSeries, RasterSeries
export Projected, Mapped
export Band
export missingval, boolmask, missingmask, replace_missing, replace_missing!,
       aggregate, aggregate!, disaggregate, disaggregate!, mask, mask!,
       resample, warp, zonal, crop, extend, trim, slice, combine, points,
       classify, classify!, mosaic, mosaic!, extract, rasterize, rasterize!,
       coverage, coverage!, setcrs, setmappedcrs, smapseries, cellsize
export crs, mappedcrs, mappedindex, mappedbounds, projectedindex, projectedbounds
export reproject, convertlookup


const DD = DimensionalData
const DA = DiskArrays
const GI = GeoInterface
const LA = Lookups

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

const DEFAULT_POINT_ORDER = (X(), Y())
const DEFAULT_TABLE_DIM_KEYS = (:X, :Y)

include("nokw.jl")
include("methods/shared_docstrings.jl")
include("lookup.jl")
include("dimensions.jl")
include("sources/sources.jl")
include("modifieddiskarray.jl")
include("filearray.jl")
include("filestack.jl")
include("openstack.jl")
include("array.jl")
include("stack.jl")
include("series.jl")
include("crs.jl")

const RasterStackOrArray = Union{AbstractRasterStack,AbstractRaster}
const RasterSeriesOrStack = Union{AbstractRasterSeries,AbstractRasterStack}

include("utils.jl")
include("skipmissing.jl")

include("table_ops.jl")
include("create.jl")
include("read.jl")
include("write.jl")
include("show.jl")
include("plotrecipes.jl")

include("methods/burning/edges.jl")
include("methods/burning/allocs.jl")
include("methods/burning/array_init.jl")
include("methods/burning/geometry.jl")
include("methods/burning/point.jl")
include("methods/burning/line.jl")
include("methods/burning/polygon.jl")
include("methods/burning/extents.jl")
include("methods/burning/utils.jl")

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
include("methods/slice_combine.jl")
include("methods/trim.jl")
include("methods/zonal.jl")

include("sources/grd.jl")
include("sources/commondatamodel.jl")
include("extensions.jl")

# Compatibility with pre-1.9 julia
function __init__()
    @static if !isdefined(Base, :get_extension)
        @require ArchGDAL = "c9ce4bd3-c3d5-55b8-8973-c0e20141b8c3" include("../ext/RastersArchGDALExt/RastersArchGDALExt.jl")
        @require CoordinateTransformations = "150eb455-5306-5404-9cee-2592286d6298" include("../ext/RastersCoordinateTransformationsExt/RastersCoordinateTransformationsExt.jl")
        @require HDF5 = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f" include("../ext/RastersHDF5Ext/RastersHDF5Ext.jl")
        @require Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("../ext/RastersMakieExt/RastersMakieExt.jl")
        @require NCDatasets = "85f8d34a-cbdd-5861-8df4-14fed0d494ab" include("../ext/RastersNCDatasetsExt/RastersNCDatasetsExt.jl")
        @require GRIBDatasets = "82be9cdb-ee19-4151-bdb3-b400788d9abc" include("../ext/RastersGRIBDatasetsExt/RastersGRIBDatasetsExt.jl")
        @require RasterDataSources = "3cb90ccd-e1b6-4867-9617-4276c8b2ca36" include("../ext/RastersRasterDataSourcesExt/RastersRasterDataSourcesExt.jl")
    end
end

end
