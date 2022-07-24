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

# Example

```julia
using Rasters, Shapefile, Downloads, Statistics, DataFrames

# Download a borders shapefile
shp_url = "https://github.com/nvkelso/natural-earth-vector/raw/master/10m_cultural/ne_10m_admin_0_countries.shp"
dbf_url = "https://github.com/nvkelso/natural-earth-vector/raw/master/10m_cultural/ne_10m_admin_0_countries.dbf"
shp_name, dbf_name = "country_borders.shp", "country_borders.dbf"
isfile(shapefile_name) || Downloads.download(shp_url, shp_name)
isfile(dbf_url) || Downloads.download(dbf_url, dbf_name)

# Download and read a raster stack
st = RasterStack(WorldClim{Climate}; month=Jan, lazy=false)

# Load the shapes for world countries
countries = Shapefile.Table(shapefile_name) |> DataFrame
# Calculate the january mean of all climate variables for all countries
january_stats = zonal(mean, st, countries) |> DataFrame
# Add the country name column (has some errors in the dataset strings it seems)
insertcols!(january_stats, 1, :country => countries.ADMIN)
```
"""
zonal(f, x::RasterStackOrArray; of) = _zonal(f, x, of)

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
_zonal(f, x::RasterStackOrArray, of) = _zonal(f, x, GI.trait(of), of)
function _zonal(f, x, ::GI.AbstractFeatureCollectionTrait, fc)
    [_zonal(f, x, feature) for feature in GI.getfeature(fc)]
end
_zonal(f, x::RasterStackOrArray, ::GI.AbstractFeatureTrait, feature) = _zonal(f, x, GI.geometry(feature))
function _zonal(f, x::AbstractRaster, ::GI.AbstractGeometryTrait, geom)
    zone = mask(crop(x; to=geom); with=geom)
    prod(size(zone)) > 0 || return missing
    return f(skipmissing(zone))
end
function _zonal(f, x::AbstractRasterStack, ::GI.AbstractGeometryTrait, geom)
    zone = mask(crop(x; to=geom); with=geom)
    return map(zone) do A
        prod(size(A)) > 0 || return missing
        f(skipmissing(A))
    end
end
function _zonal(f, x::RasterStackOrArray, ::Nothing, obj)
    if Tables.istable(obj)
        geoms = Tables.getcolumn(obj, first(GI.geometrycolumns(obj)))
        return [_zonal(f, x, geom) for geom in geoms]
    elseif obj isa AbstractVector
        return [_zonal(f, x, geom) for geom in obj]
    else
        throw(ArgumentError("Cannot calculate zonal statistics for objects of type $(typeof(obj))"))
    end
end
