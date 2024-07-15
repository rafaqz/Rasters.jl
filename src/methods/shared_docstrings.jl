# Share common docstrings here to keep things consistent

const TO_KEYWORD = """
- `to`: a `Raster`, `RasterStack`, `Tuple` of `Dimension` or `Extents.Extent`.
    If no `to` object is provided the extent will be calculated from the geometries,
    Additionally, when no `to` object or an `Extent` is passed for `to`, the `size`
    or `res` keyword must also be used.
"""
const SIZE_KEYWORD = """
- `size`: the size of the output array, as a `Tuple{Int,Int}` or single `Int` for a square.
    Only required when `to` is not used or is an `Extents.Extent`, and `res` is not used.
"""
const RES_KEYWORD = """
- `res`: the resolution of the dimensions (often in meters or degrees), a `Real` or `Tuple{<:Real,<:Real}`.
    Only required when `to` is not used or is an `Extents.Extent`, and `size` is not used.
"""
const CRS_KEYWORD = """
- `crs`: a `crs` which will be attached to the resulting raster when `to` not passed
   or is an `Extent`. Otherwise the crs from `to` is used.
"""

const SHAPE_KEYWORDS = """
- `shape`: Force `data` to be treated as `:polygon`, `:line` or `:point` geometries.
    using points or lines as polygons may have unexpected results.
- `boundary`: for polygons, include pixels where the `:center` is inside the polygon,
    where the polygon `:touches` the pixel, or that are completely `:inside` the polygon.
    The default is `:center`.
"""

const THREADED_KEYWORD = """
- `threaded`: run operations in parallel, `false` by default. In some circumstances `threaded` 
    can give large speedups over single-threaded operation. This can be true for complicated 
    geometries written into low-resolution rasters, but may not be for simple geometries with 
    high-resolution rasters. With very large rasters threading may be counter productive due 
    to excessive memory use. Caution should also be used: `threaded` should not be used in in-place 
    functions writing to `BitArray` or other arrays where race conditions can occur. 
"""

const GEOM_KEYWORDS = """
$TO_KEYWORD
$RES_KEYWORD
$SIZE_KEYWORD
$CRS_KEYWORD
$SHAPE_KEYWORDS
"""

const DATA_ARGUMENT = """
- `data`: a GeoInterface.jl `AbstractGeometry`, a nested `Vector` of `AbstractGeometry`,
    or a Tables.jl compatible object containing a `:geometry` column or points and values columns,
    in which case `geometrycolumn` must be specified.
"""

const OBJ_ARGUMENT = """
a [`Raster`](@ref) or one or multiple geometries. Geometries can be
    a GeoInterface.jl `AbstractGeometry`, a nested `Vector` of `AbstractGeometry`,
    or a Tables.jl compatible object containing a `:geometry` column or points and values columns,
    in which case `geometrycolumn` must be specified."""

const GEOMETRYCOLUMN_KEYWORD = """
- `geometrycolumn`: `Symbol` to manually select the column the geometries are in
    when `data` is a Tables.jl compatible table, or a tuple of `Symbol` for columns of
    point coordinates.
"""


const LAZY_KEYWORD = """
- `lazy`: A `Bool` specifying if to load data lazily from disk. `false` by default.
"""
const FILENAME_KEYWORD = """
- `filename`: a filename to write to directly, useful for large files.
"""
const SUFFIX_KEYWORD = """
- `suffix`: a string or value to append to the filename.
    A tuple of `suffix` will be applied to stack layers. `keys(stack)` are the default.
"""
const PROGRESS_KEYWORD = """
- `progress`: show a progress bar, `true` by default, `false` to hide.
"""
const VERBOSE_KEYWORD = """
- `vebose`: whether to print messages about potential problems. `true` by default.
"""
const SOURCE_KEYWORD = """
- `source`: Usually automatically detected from filepath extension. 
    To manually force, a `Symbol` can be passed `:gdal`, `:netcdf`, `:grd`, `:grib`.
    The internal [`Rasters.Source`](@ref) objects, such as `Rasters.GDALsource()`, 
    `Rasters.GRIBsource()` or `Rasters.NCDsource()` can also be used.
"""
const EXT_KEYWORD = """
- `ext`: filename extension such as ".tiff" or ".nc". 
    Used to specify specific files if only a directory path is used.
"""
const FORCE_KEYWORD = """
- `force`: `false` by default. If `true` it force writing to a file destructively, even if it already exists.
"""
const DROPBAND_KEYWORD = """
- `dropband`: drop single band dimensions when creating stacks from filenames. `true` by default.
"""

const CONSTRUCTOR_CRS_KEYWORD = """
- `crs`: the coordinate reference system of  the objects `XDim`/`YDim` dimensions.
    Only set this if you know the detected crs is incorrect, or it is not present in
    the file. The `crs` is expected to be a GeoFormatTypes.jl `CRS` or `Mixed` mode `GeoFormat` object,
    like `EPSG(4326)`.
"""
const CONSTRUCTOR_MAPPEDCRS_KEYWORD = """
- `mappedcrs`: the mapped coordinate reference system of the objects `XDim`/`YDim` dimensions.
    for `Mapped` lookups these are the actual values of the index. For `Projected` lookups
    this can be used to index in eg. `EPSG(4326)` lat/lon values, having it converted automatically.
    Only set this if the detected `mappedcrs` in incorrect, or the file does not have a `mappedcrs`,
    e.g. a tiff. The `mappedcrs` is expected to be a GeoFormatTypes.jl `CRS` or `Mixed` mode `GeoFormat` type.
"""
const GROUP_KEYWORD = """
- `group`: the group in the dataset where `name` can be found. Only needed for nested datasets.
    A `String` or `Symbol` will select a single group. Pairs can also used to access groups
    at any nested depth, i.e `group=:group1 => :group2 => :group3`.
"""

const REPLACE_MISSING_KEYWORD = """
- `replace_missing`: replace `missingval` with `missing`. This is done lazily if `lazy=true`.
    Note that currently for NetCDF and GRIB files `replace_missing` is always true. 
    In future `replace_missing=false` will also work for these data sources.
"""

const CHECKMEMORY_KEYWORD = """
- `checkmemory`: If `true` (the default), check if there is enough memory for the operation. 
    `false` will ignore memory needs.
"""
