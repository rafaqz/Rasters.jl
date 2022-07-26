"""
    zonal(f, x::RasterStackOrArray; of)

Calculate zonal statistics for the the zone of a `Raster` or `RasterStack` defined
the `of` keyword, which can be an `Extent`, or GeoInterface.jl compatible object. 
Tables, `AbstractVectors` an feature collections will return `Vector`s of values.

# Arguments

- `f`: any function that reduces an iterable to a single value, such as `Base.sum` or `Statistic.mean`
- `x`: A `Raster` or `RasterStack`
- `zone`: A `Raster`, `RasterStack`, dim tuple, extent, GeoInterface.jl compatible geometry,
    Tables.jl compatible table of a `:geometry` column, or an `AbstractVector` of
    any of these objects..

# Keywords

These can be used when `of` is a GeoInterface.jl compatible object:

- `shape`: Force `data` to be treated as `:polygon`, `:line` or `:point`, where possible.
- `boundary`: for polygons, include pixels where the `:center` is inside the polygon,
    where the line `:touches` the pixel, or that are completely `:inside` inside the polygon.
    The default is `:center`.

# Example

```jldoctest
using Rasters, Shapefile, DataFrames, Downloads, Statistics, Dates

# Download a borders shapefile
ne_url = "https://github.com/nvkelso/natural-earth-vector/raw/master/10m_cultural/ne_10m_admin_0_countries"
shp_url, dbf_url  = ne_url * ".shp", ne_url * ".dbf"
shp_name, dbf_name = "country_borders.shp", "country_borders.dbf"
isfile(shp_name) || Downloads.download(shp_url, shp_name)
isfile(dbf_url) || Downloads.download(dbf_url, dbf_name)

# Download and read a raster stack from WorldClim
st = RasterStack(WorldClim{Climate}; month=Jan, lazy=false)

# Load the shapes for world countries
countries = Shapefile.Table(shp_name) |> DataFrame
# Calculate the january mean of all climate variables for all countries
january_stats = zonal(mean, st; of=countries, boundary=:touches) |> DataFrame
# Add the country name column (natural earth has some string errors it seems)
insertcols!(january_stats, 1, :country => countries.ADMIN)

# output
"""
zonal(f, x::RasterStackOrArray; of, shape=nothing, boundary=nothing) = 
    _zonal(f, x, of; shape, boundary)

_zonal(f, x::RasterStackOrArray, of::RasterStackOrArray) = _zonal(f, x, Extents.extent(of))
_zonal(f, x::RasterStackOrArray, of::DimTuple) = _zonal(f, x, Extents.extent(of))
# We don't need to `mask` with an extent, it's square so `crop` will do enough.
_zonal(f, x::Raster, of::Extents.Extent) = f(skipmissing(crop(x; to=of)))
function _zonal(f, x::RasterStack, ext::Extents.Extent)
    zone = crop(x; to=ext)
    prod(size(zone)) > 0 || return missing
    return map(zone) do A
        f(skipmissing(A))
    end
end
# Otherwise of is a geom, table or vector
_zonal(f, x::RasterStackOrArray, of; kw...) = _zonal(f, x, GI.trait(of), of; kw...)
function _zonal(f, x, ::GI.AbstractFeatureCollectionTrait, fc; kw...)
    [_zonal(f, x, feature; kw...) for feature in GI.getfeature(fc)]
end
_zonal(f, x::RasterStackOrArray, ::GI.AbstractFeatureTrait, feature; kw...) =
    _zonal(f, x, GI.geometry(feature); kw...)
function _zonal(f, x::AbstractRaster, ::GI.AbstractGeometryTrait, geom; kw...)
    cropped = crop(x; to=geom)
    prod(size(cropped)) > 0 || return missing
    zone = mask(cropped; with=geom, kw...)
    return f(skipmissing(zone))
end
function _zonal(f, st::AbstractRasterStack, ::GI.AbstractGeometryTrait, geom; kw...)
    cropped = crop(st; to=geom)
    prod(size(first(cropped))) > 0 || return map(_ -> missing, st)
    masked = mask(cropped; with=geom, kw...)
    return map(masked) do A
        prod(size(A)) > 0 || return missing
        f(skipmissing(A))
    end
end
function _zonal(f, x::RasterStackOrArray, ::Nothing, obj; kw...)
    if Tables.istable(obj)
        geoms = Tables.getcolumn(obj, first(GI.geometrycolumns(obj)))
        return [_zonal(f, x, geom; kw...) for geom in geoms]
    elseif obj isa AbstractVector
        return [_zonal(f, x, geom; kw...) for geom in obj]
    else
        throw(ArgumentError("Cannot calculate zonal statistics for objects of type $(typeof(obj))"))
    end
end
