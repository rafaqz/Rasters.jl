module GeoData

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end GeoData

using Dates

import Adapt,
       ArchGDAL,
       ColorTypes,
       ConstructionBase,
       DiskArrays,
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
       Setfield

Reexport.@reexport using DimensionalData, GeoFormatTypes, RasterDataSources

using DimensionalData.Tables

using RecipesBase: @recipe, @series
using Base: tail, @propagate_inbounds
using DimensionalData: StandardIndices, DimTuple
using Setfield: @set, @set!
using ColorTypes: RGB

export AbstractGeoArray, GeoArray
export AbstractGeoStack, GeoStack
export AbstractGeoSeries, GeoSeries
export Projected, Mapped
export Band
export missingval, boolmask, missingmask, replace_missing,
       aggregate, aggregate!, disaggregate, disaggregate!, mask, mask!, 
       resample, warp, crop, extend, trim, slice, points, subset, inpolygon,
       classify, classify!, mosaic, mosaic!, extract, rasterize, rasterize!
export crs, mappedcrs, mappedindex, mappedbounds, projectedindex, projectedbounds
export reproject, convertmode
export geoarray, stack, series


const DD = DimensionalData
const DA = DiskArrays
const GI = GeoInterface

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

include("mode.jl")
include("dimensions.jl")
include("filearray.jl")
include("array.jl")
include("filestack.jl")
include("stack.jl")
include("series.jl")
include("utils.jl")
include("aggregate.jl")
include("methods.jl")
include("read.jl")
include("sources/grd.jl")
include("show.jl")
include("plotrecipes.jl")
include("convenience.jl")
include("write.jl")

include("sources/smap.jl")
include("sources/ncdatasets.jl")
include("sources/gdal.jl")
include("reproject.jl")
include("sources/rasterdatasources.jl")

end
