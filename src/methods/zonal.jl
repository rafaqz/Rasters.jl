"""
    zonal(f, x::Union{Raster,RasterStack}; of, kw...)

Calculate zonal statistics for the the zone of a `Raster` or `RasterStack`
covered by the `of` object/s.

# Arguments

- `f`: any function that reduces an iterable to a single value, such as `sum` or `Statistics.mean`
- `x`: A `Raster` or `RasterStack`
- `of`: A `DimTuple`, `Extent`, $OBJ_ARGUMENT

# Keywords
$GEOMETRYCOLUMN_KEYWORD
These can be used when `of` is or contains (a) GeoInterface.jl compatible object(s):

- `shape`: Force `data` to be treated as `:polygon`, `:line` or `:point`, where possible.
- `boundary`: for polygons, include pixels where the `:center` is inside the polygon,
    where the line `:touches` the pixel, or that are completely `:inside` inside the polygon.
    The default is `:center`.
- `progress`: show a progress bar, `true` by default, `false` to hide..
- `skipmissing`: wether to apply `f` to the result of `skipmissing(A)` or not. If `true`
    `f` will be passed an iterator over the values, which loses all spatial information.
    if `false` `f` will be passes a masked `Raster` or `RasterStack`, and will be responsible
    for handling missing values itself. The default value is `true`.

# Example

```jldoctest
using Rasters, RasterDataSources, ArchGDAL, Shapefile, DataFrames, Downloads, Statistics, Dates

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
january_stats = zonal(mean, st; of=countries, boundary=:touches, progress=false) |> DataFrame

# Add the country name column (natural earth has some string errors it seems)
insertcols!(january_stats, 1, :country => first.(split.(countries.ADMIN, r"[^A-Za-z ]")))

# output
258×8 DataFrame
 Row │ country                       tmin       tmax       tavg       prec     ⋯
     │ SubStrin…                     Float32    Float32    Float32    Float64  ⋯
─────┼──────────────────────────────────────────────────────────────────────────
   1 │ Indonesia                      21.5447    29.1864    25.3656   271.063  ⋯
   2 │ Malaysia                       21.3087    28.4291    24.8688   273.381
   3 │ Chile                           7.24534   17.9263    12.5858    78.1287
   4 │ Bolivia                        17.2065    27.7454    22.4759   192.542
   5 │ Peru                           15.0273    25.5504    20.2888   180.007  ⋯
   6 │ Argentina                      13.6751    27.6715    20.6732    67.1837
   7 │ Dhekelia Sovereign Base Area    5.87126   15.8991    10.8868    76.25
   8 │ Cyprus                          5.65921   14.6665    10.1622    97.4474
  ⋮  │              ⋮                    ⋮          ⋮          ⋮         ⋮     ⋱
 252 │ Spratly Islands                25.0       29.2       27.05      70.5    ⋯
 253 │ Clipperton Island              21.5       33.2727    27.4        6.0
 254 │ Macao S                        11.6694    17.7288    14.6988    28.0
 255 │ Ashmore and Cartier Islands   NaN        NaN        NaN        NaN
 256 │ Bajo Nuevo Bank               NaN        NaN        NaN        NaN      ⋯
 257 │ Serranilla Bank               NaN        NaN        NaN        NaN
 258 │ Scarborough Reef              NaN        NaN        NaN        NaN
                                                  3 columns and 243 rows omitted
```
"""
zonal(f, x::RasterStackOrArray; of, kw...) = _zonal(f, x, of; kw...)

_zonal(f, x::RasterStackOrArray, of::RasterStackOrArray; kw...) = 
    _zonal(f, x, Extents.extent(of); kw...)
_zonal(f, x::RasterStackOrArray, of::DimTuple; kw...) = 
    _zonal(f, x, Extents.extent(of); kw...)
# We don't need to `mask` with an extent, it's square so `crop` will do enough.
_zonal(f, x::Raster, of::Extents.Extent; skipmissing=true) =
    _maybe_skipmissing_call(f, crop(x; to=of, touches=true), skipmissing)
function _zonal(f, x::RasterStack, ext::Extents.Extent; skipmissing=true)
    cropped = crop(x; to=ext, touches=true)
    prod(size(cropped)) > 0 || return missing
    return maplayers(cropped) do A
        _maybe_skipmissing_call(f, A, skipmissing)
    end
end
# Otherwise of is a geom, table or vector
_zonal(f, x::RasterStackOrArray, of; kw...) = _zonal(f, x, GI.trait(of), of; kw...)

_zonal(f, x, ::GI.AbstractFeatureCollectionTrait, fc; kw...) =
    _zonal(f, x, nothing, fc; kw...)
_zonal(f, x::RasterStackOrArray, ::GI.AbstractFeatureTrait, feature; kw...) =
    _zonal(f, x, GI.geometry(feature); kw...)
function _zonal(f, x::AbstractRaster, ::GI.AbstractGeometryTrait, geom; 
    skipmissing=true, kw...
)
    cropped = crop(x; to=geom, touches=true)
    prod(size(cropped)) > 0 || return missing
    # Fast path for lines
    if skipmissing && trait isa GI.AbstractLineStringTrait
        # Read the lines directly and calculate stats on the fly
        return _zonal_lines(f, _reduce_op(f), st, geom; kw...)
    else
        # Mask the cropped raster 
        # TODO use mask! on re-used buffer where possible
        masked = mask(cropped; with=geom, kw...)
        # Calculate stats on whole raster or `skipmissing(raster)``    
        return _maybe_skipmissing_call(f, masked, skipmissing)
    end
end
function _zonal(f, st::AbstractRasterStack, trait::GI.AbstractGeometryTrait, geom; 
    skipmissing=true, kw...
)
    cropped = crop(st; to=geom, touches=true)
    prod(size(cropped)) > 0 || return map(_ -> missing, layerdims(st))
    # Fast path for lines
    if skipmissing && trait isa GI.AbstractLineStringTrait
        # Read the lines directly and calculate stats on the fly
        return _zonal_lines(f, _reduce_op(f), cropped, geom; kw...)
    else
        # Mask the cropped stack 
        # TODO use mask! on re-used buffer where possible
        masked = mask(cropped; with=geom, kw...)
        # Calculate stats on whole stack or `skipmissing(raster)``    
        return maplayers(masked) do A
            prod(size(A)) > 0 || return missing
            _maybe_skipmissing_call(f, A, skipmissing)
        end
    end
end
function _zonal(f, x::RasterStackOrArray, ::Nothing, data; 
    geometrycolumn=nothing, kw...
)
    geoms = _get_geometries(data, geometrycolumn)
    _zonal(f, x, nothing, geoms; kw...) 
end
function _zonal(f, x::RasterStackOrArray, ::Nothing, geoms::AbstractArray; 
    progress=true, threaded=true, geometrycolumn=nothing, kw...
)
    n = length(geoms)
    n == 0 && return [] # TODO make this type stable
    # Allocate vector by running _zonal on the first geoms
    # until a non-missing value is returns
    zs, start_index = _alloc_zonal(f, x, geoms, n; kw...)
    # If all geometries are missing, return a vector of missings
    start_index == n + 1 && return zs
    # Otherwise run _zonal for the rest of the geometries, possibly threaded
    _run(start_index:n, threaded, progress, "Applying $f to each geometry...") do i
        zs[i] = _zonal(f, x, geoms[i]; kw...)
    end
    return zs
end

function _zonal_lines(f, op, r::RasterStackOrArray, geom; 
    buffer, kw...
)
    burn_lines(rast, geom) do i
        buffer[i] = rast[i]
    end
end
function _zonal_lines(f, op::Nothing, st::RasterStackOrArray, geom; kw...)
    # burn and count
end

function _alloc_zonal(f, x, geoms, n; kw...)
    # Find first non-missing entry and count number of missing entries
    n_missing::Int = 0
    z1 = missing
    # If the first stat is missing, repeat until we find one that is non-missing
    for geom in geoms
        z1 = _zonal(f, x, geom; kw...)
        ismissing(z1) || break
        n_missing += 1
    end
    zs = Vector{Union{Missing,typeof(z1)}}(undef, n)
    zs[1:n_missing] .= missing
    # Exit early when all elements are missing
    if n_missing == n
        return zs, n_missing + 1
    end
    zs[n_missing + 1] = z1
    return zs, n_missing + 1
end

_maybe_skipmissing_call(f, A, sm) = sm ? f(skipmissing(A)) : f(A)
