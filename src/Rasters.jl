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
       GeometryOps,
       GeometryOpsCore,
       OffsetArrays,
       ProgressMeter,
       Missings,
       Mmap,
       RecipesBase,
       Reexport,
       Setfield,
       SortTileRecursiveTree,
       Statistics

Reexport.@reexport using DimensionalData, GeoFormatTypes

using DimensionalData.Tables,
      DimensionalData.Lookups,
      DimensionalData.Dimensions,
      DimensionalData.Lookups.IntervalSets

using DimensionalData: Name, NoName
using .Dimensions: StandardIndices, DimTuple
using .Lookups: LookupTuple

using Statistics: mean
using RecipesBase: @recipe, @series
using Base: tail, @propagate_inbounds

import GeoInterface: crs
import Extents: Extent, extent

using Setfield: @set, @set!
using ColorTypes: RGB

using CommonDataModel: AbstractDataset, AbstractVariable

using DiskArrays: @implement_diskarray, eachchunk, haschunks, isdisk

using GeometryOpsCore: Planar, Spherical
export Planar, Spherical

export AbstractRaster, Raster
export AbstractRasterStack, RasterStack
export AbstractRasterSeries, RasterSeries
export Projected, Mapped, GeometryLookup
export Band, Geometry
export missingval, boolmask, missingmask, replace_missing, replace_missing!,
       aggregate, aggregate!, disaggregate, disaggregate!, mask, mask!,
       resample, warp, zonal, crop, extend, trim, slice, combine, points,
       classify, classify!, mosaic, mosaic!, extract, rasterize, rasterize!,
       coverage, coverage!, setcrs, setmappedcrs, smapseries, cellsize, cellarea
export crs, mappedcrs, mappedindex, mappedbounds, projectedindex, projectedbounds
export reproject, convertlookup
export Extent, extent

const DD = DimensionalData
const DA = DiskArrays
const GI = GeoInterface
const GO = GeometryOps
const LA = Lookups

# DimensionalData documentation urls
const DDdocs = "https://rafaqz.github.io/DimensionalData.jl/stable/api/reference"
const DDdimdocs = join([DDdocs, "#DimensionalData.Dimensions.Dimension"])
const DDarraydocs = join([DDdocs, "#DimensionalData.AbstractDimArray"])
const DDabssampleddocs = join([DDdocs, "#DimensionalData.Dimensions.Lookups.AbstractSampled"])
const DDsampleddocs = join([DDdocs, "#DimensionalData.Dimensions.Lookups.Sampled"])
const DDlocusdocs = join([DDdocs, "#DimensionalData.Dimensions.Lookups.locus"])
const DDselectordocs = join([DDdocs, "#DimensionalData.Dimensions.Lookups.Selector"])
const DDtidocs = join([DDdocs, "#DimensionalData.Ti"])
const DDregulardocs = join([DDdocs, "#DimensionalData.Lookups.Regular"])
const DDirregulardocs = join([DDdocs, "#DimensionalData.Lookups.Irregular"])


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

include("geometry_lookup/geometry_lookup.jl")
include("geometry_lookup/lookups.jl")
include("geometry_lookup/methods.jl")
include("geometry_lookup/io.jl")

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

end
