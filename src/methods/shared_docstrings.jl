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
    to excessing memory use. Caution should also be used: `threaded` should not be used in in-place 
    functions wrinting to `BitArray` or other arrays where race conditions can occur. 
"""

const GEOM_KEYWORDS = """
$TO_KEYWORD
$RES_KEYWORD
$SIZE_KEYWORD
$CRS_KEYWORD
$SHAPE_KEYWORDS
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
