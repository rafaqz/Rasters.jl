"""
    zonal(f, x::RasterStackOrArray, zone)

Calculate zonal statistics for the the zone of a `Raster` or `RasterStack` defined
by the object `zone`.

# Arguments

- `f`: aby function that reduces an iterable to a single value, such as `Base.sum` or `Statistic.mean`
- `x`: A `Raster` or `RasterStack`
- `zone`: A `Raster`, `RasterStack`, dim tuple, extent, GeoInterface.jl compatible geometry,
    Tables.jl compatible table with a `:geometry` column, or an `AbstractVector` of
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
function zonal end
zonal(f, x::RasterStackOrArray, obj::RasterStackOrArray) = zonal(f, x, Extents.extent(obj))
zonal(f, x::RasterStackOrArray, dims::DimTuple) = zonal(f, x, Extents.extent(dims))
zonal(f, x::Raster, obj::Extents.Extent) = f(skipmissing(crop(x; to=obj)))
function zonal(f, x::RasterStack, ext::Extents.Extent) 
    zone = crop(x; to=ext)
    prod(size(zone)) > 0 || return missing
    return map(zone) do A
        f(skipmissing(A))
    end
end
zonal(f, x::RasterStackOrArray, obj) = _zonal(f, x, GI.trait(obj), obj)

function _zonal(f, x, ::GI.AbstractFeatureCollectionTrait, fc)
    [zonal(f, x, feature) for feature in GI.getfeature(fc)]
end
_zonal(f, x, ::GI.AbstractFeatureTrait, obj) = _zonal(f, x, GI.geometry(obj))
function _zonal(f, x::AbstractRaster, ::GI.AbstractGeometryTrait, obj)
    zone = crop(x; to=obj)
    prod(size(zone)) > 0 || return missing
    return f(skipmissing(zone))
end
function _zonal(f, x::AbstractRasterStack, ::GI.AbstractGeometryTrait, obj)
    zone = crop(x; to=obj)
    return map(zone) do A
        prod(size(A)) > 0 || return missing
        f(skipmissing(A))
    end
end
function _zonal(f, x, ::Nothing, obj)
    if Tables.istable(obj)
        geoms = Tables.getcolumn(obj, first(GI.geometrycolumns(obj)))
        return [zonal(f, x, geom) for geom in geoms]
    elseif obj isa AbstractVector
        return [zonal(f, x, geom) for geom in obj]
    else
        throw(ArgumentError("Cannot calculate zonal statistics for objects of type $(typeof(obj))"))
    end
end
