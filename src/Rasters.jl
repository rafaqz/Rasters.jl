module Rasters

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end Rasters

using Dates

import Adapt,
       ArchGDAL,
       ColorTypes,
       CoordinateTransformations,
       ConstructionBase,
       DiskArrays,
       Extents,
       FillArrays,
       Flatten,
       GeoInterface,
       HDF5,
       PolygonInbounds,
       ProgressMeter,
       Missings,
       Mmap,
       NCDatasets,
       RecipesBase,
       Reexport,
       Setfield,
       ASCIIrasters

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
       resample, warp, zonal, crop, extend, trim, slice, points, subset, inpolygon,
       classify, classify!, mosaic, mosaic!, extract, rasterize, rasterize!,
       setcrs, setmappedcrs
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

# Source dispatch singletons
struct NCDfile end
struct GRDfile end
struct GDALfile end
struct SMAPfile end
struct ASCIIfile end

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
include("convenience.jl")
include("show.jl")
include("plotrecipes.jl")


include("methods/aggregate.jl")
include("methods/classify.jl")
include("methods/crop_extend.jl")
include("methods/extract.jl")
include("methods/inpolygon.jl")
include("methods/mask.jl")
include("methods/mosaic.jl")
include("methods/points.jl")
include("methods/rasterize.jl")
include("methods/replace_missing.jl")
include("methods/reproject.jl")
include("methods/resample.jl")
include("methods/slice_combine.jl")
include("methods/trim.jl")
include("methods/warp.jl")
include("methods/zonal.jl")

include("sources/grd.jl")
include("sources/smap.jl")
include("sources/ncdatasets.jl")
include("sources/ascii.jl")
include("sources/gdal.jl")
include("sources/rasterdatasources.jl")

end
