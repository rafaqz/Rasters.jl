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
- `res`: the resolution of the dimensions, a `Real` or `Tuple{<:Real,<:Real}`.
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

const GEOM_KEYWORDS = """
$TO_KEYWORD
$RES_KEYWORD
$SIZE_KEYWORD
$CRS_KEYWORD
$SHAPE_KEYWORDS
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
