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
- `spatialslices`: if `true`, and the input Raster has more dimensions than `X` and `Y`, then we will 
    apply the function `f` to each spatial slice of the raster (literally, `mapslices(f, x; dims = (X, Y))`),
    and return a vector data cube or stack of the results.
    if `false`, then we will apply `f` to the full cropped raster, and return a vector of results (one per geometry)
    as usual.

# Example

```jldoctest
using Rasters, RasterDataSources, ArchGDAL, DataFrames, Statistics, Dates, NaturalEarth
# Download borders
countries = naturalearth("admin_0_countries", 10) |> DataFrame
# Download and read a raster stack from WorldClim
st = RasterStack(WorldClim{Climate}; month=Jan)
# Calculate the january mean of all climate variables for all countries
january_stats = zonal(mean, st; of=countries, boundary=:touches, progress=false) |> DataFrame
# Add the country name column (natural earth has some string errors it seems)
insertcols!(january_stats, 1, :country => first.(split.(countries.ADMIN, r"[^A-Za-z ]")))
# output
258×8 DataFrame
 Row │ country                       tmin       tmax       tavg       prec     ⋯
     │ SubStrin…                     Float32    Float32    Float32    Float64  ⋯
─────┼──────────────────────────────────────────────────────────────────────────
   1 │ Indonesia                      21.5447    29.1865    25.3656   271.063  ⋯
   2 │ Malaysia                       21.3087    28.4291    24.8688   273.381
   3 │ Chile                           7.24534   17.9262    12.5858    78.1287
   4 │ Bolivia                        17.2065    27.7454    22.4758   192.542
   5 │ Peru                           15.0273    25.5504    20.2889   180.007  ⋯
   6 │ Argentina                      13.6751    27.6716    20.6732    67.1837
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
function zonal(f, x::RasterStack; of, skipmissing=true, spatialslices=_True(), kw...)
    # TODO: open currently doesn't work so well for large rasterstacks,
    # we need to fix that before we can go back to this being a single method
    # on `RasterStackOrArray`.
    _zonal(f, _prepare_for_burning(x), of; skipmissing, spatialslices, kw...)
end
function zonal(f, x::Raster; of, skipmissing=true, spatialslices=_True(), kw...)
    open(x) do xo
        _zonal(f, _prepare_for_burning(xo), of; skipmissing, spatialslices, kw...)
    end
end

_zonal(f, x::RasterStackOrArray, of::RasterStackOrArray; kw...) = 
    _zonal(f, x, Extents.extent(of); kw...)
_zonal(f, x::RasterStackOrArray, of::DimTuple; kw...) = 
    _zonal(f, x, Extents.extent(of); kw...)
# We don't need to `mask` with an extent, it's square so `crop` will do enough.
_zonal(f, x::Raster, of::Extents.Extent; skipmissing, spatialslices) = _maybe_skipmissing_call(_maybe_spatialsliceify(f, spatialslices), crop(x; to=of, touches=true), skipmissing)
function _zonal(f, x::RasterStack, ext::Extents.Extent; skipmissing, spatialslices)
    cropped = crop(x; to=ext, touches=true)
    if length(cropped) == 0 && skipmissing == true
        return map(_ -> missing, x)
    end
    return maplayers(cropped) do A
        _maybe_skipmissing_call(_maybe_spatialsliceify(f, spatialslices), A, skipmissing)
    end
end
# Otherwise of is a geom, table or vector
_zonal(f, x::RasterStackOrArray, of; kw...) = _zonal(f, x, GI.trait(of), of; kw...)

_zonal(f, x, ::GI.AbstractFeatureCollectionTrait, fc; kw...) =
    _zonal(f, x, nothing, fc; kw...)
_zonal(f, x::RasterStackOrArray, ::GI.AbstractFeatureTrait, feature; kw...) =
    _zonal(f, x, GI.geometry(feature); kw...)
function _zonal(f, x::AbstractRaster, ::GI.AbstractGeometryTrait, geom; 
    skipmissing, spatialslices, kw...
)
    cropped = crop(x; to=geom, touches=true)
    if length(cropped) == 0 && skipmissing == true
        return missing
    end
    masked = mask(cropped; with=geom, kw...)
    return _maybe_skipmissing_call(_maybe_spatialsliceify(f, spatialslices), masked, skipmissing)
end
function _zonal(f, st::AbstractRasterStack, ::GI.AbstractGeometryTrait, geom; 
    skipmissing, spatialslices, kw...
)
    cropped = crop(st; to=geom, touches=true)
    if length(cropped) == 0 && skipmissing == true
        return map(_ -> missing, st)
    end
    masked = mask(cropped; with=geom, kw...)
    return maplayers(masked) do A
        if length(A) == 0 && skipmissing == true
            return missing
        end
        _maybe_skipmissing_call(_maybe_spatialsliceify(f, spatialslices), A, skipmissing)
    end
end
function _zonal(f, x::RasterStackOrArray, ::Nothing, data; 
    progress=true, threaded=true, geometrycolumn=nothing, kw...
)
    geoms = _get_geometries(data, geometrycolumn)
    n = length(geoms)
    n == 0 && return []
    zs, start_index = _alloc_zonal(f, x, geoms, n; kw...)
    start_index == n + 1 && return zs
    _run(start_index:n, threaded, progress, "Applying $f to each geometry...") do i
        zs[i] = _zonal(f, x, geoms[i]; kw...)
    end
    return zs
end

function _alloc_zonal(f, x, geoms, n; kw...)
    # Find first non-missing entry and count number of missing entries
    n_missing::Int = 0
    z1 = _zonal(f, x, first(geoms); kw...)
    for geom in geoms
        z1 = _zonal(f, x, geom; kw...)
        if !ismissing(z1)
            break
        end
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

# Optionally wrap the input argument in `skipmissing(A)` is `sm` is true.
_maybe_skipmissing_call(f, A, sm) = istrue(sm) ? f(skipmissing(A)) : f(A)

# the only reason we have AbstractDimArray here is to make sure that DD.otherdims is available.
# We could probably get away with just AbstractArray here otherwise.
# The reason this is not just mapslices is because this drops the sliced dimensions automatically, 
# which is what we want.
function _mapspatialslices(f, x::AbstractDimArray; spatialdims = (Val{DD.XDim}(), Val{DD.YDim}()))
    iterator = DD.DimSlices(x; dims = DD.otherdims(x, spatialdims), drop = true)
    return f.(iterator)
end
# SkipMissingVal and SkipMissing both store the initial value in the `x` property,
# so we can use the same thing to extract it.
function _mapspatialslices(f, s::Union{SkipMissingVal, Base.SkipMissing}; spatialdims = (Val{DD.XDim}(), Val{DD.YDim}()))
    x = s.x # get the raster out of the SkipMissingVal
    iterator = DD.DimSlices(x; dims = DD.otherdims(x, spatialdims), drop = true)
    return @. f(skipmissing(iterator))
end

_maybe_spatialsliceify(f, spatialslices) = istrue(spatialslices) ? _SpatialSliceify(f, (Val{DD.XDim}(), Val{DD.YDim}())) : f
_maybe_spatialsliceify(f, spatialslices::DD.AllDims) = _SpatialSliceify(f, spatialslices)

"""
    _SpatialSliceify(f, dims)

A callable struct that applies `mapslices(f, x; dims = spatialdims)` to the input array `x`, and removes empty dimensions.

```jldoctest
data = ones(10, 10, 10, 10);
f = _SpatialSliceify(sum, (1, 2))
size(f(data))

# output
(10, 10)
```
"""
struct _SpatialSliceify{F, D}
    f::F
    dims::D
end

(r::_SpatialSliceify{F, D})(x) where {F, D} = _mapspatialslices(r.f, x; spatialdims = r.dims)
