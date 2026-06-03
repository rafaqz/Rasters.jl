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
- `emptyval`: value to return if a geometry returns an empty region for `f`. 
    Specifying a value for `emptyval` triggers a check `isempty` and returns `emptyval` if the result is `true`.
    If unspecified, an error will be thrown.
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
function zonal(f, x::RasterStack; of, emptyval=nokw, skipmissing=true, spatialslices=_False(), missingval=isnothing(missingval(x)) ? missing : missingval(x), kw...)
    # TODO: open currently doesn't work so well for large rasterstacks,
    # we need to fix that before we can go back to this being a single method
    # on `RasterStackOrArray`.
    _zonal(f, _prepare_for_burning(x), of; skipmissing, emptyval, spatialslices, missingval, kw...)
end
function zonal(f, x::Raster; of, emptyval=nokw, skipmissing=true, spatialslices=_False(), missingval=isnothing(missingval(x)) ? missing : missingval(x), kw...)
    open(x) do xo
        _zonal(f, _prepare_for_burning(xo), of; skipmissing, emptyval, spatialslices, missingval, kw...)
    end
end

_zonal(f, x::RasterStackOrArray, of::RasterStackOrArray; kw...) = 
    _zonal(f, x, Extents.extent(of); kw...)
_zonal(f, x::RasterStackOrArray, of::DimTuple; kw...) = 
    _zonal(f, x, Extents.extent(of); kw...)
# We don't need to `mask` with an extent, it's square so `crop` will do enough.
_zonal(f, x::Raster, of::Extents.Extent; skipmissing, emptyval, spatialslices, missingval) =
    _maybe_skipmissing_call(_maybe_spatialsliceify(f, spatialslices, missingval), crop(x; to=of, touches=true), skipmissing, emptyval)
function _zonal(f, x::RasterStack, ext::Extents.Extent; skipmissing, emptyval, spatialslices, missingval)
    cropped = crop(x; to=ext, touches=true)
    if length(cropped) == 0 && skipmissing == true
        return map(_ -> missingval, x)
    end
    return maplayers(cropped) do A
        _maybe_skipmissing_call(_maybe_spatialsliceify(f, spatialslices, missingval), A, skipmissing, emptyval)
    end
end
# GeometryLookup or dim wrapping one: forward through Geometry(...) so we go
# through the (::Nothing, ::Dimension{<:GeometryLookup}) path below.
_zonal(f, x::RasterStackOrArray, ::Nothing, data::GeometryLookup; kw...) =
    _zonal(f, x, nothing, Geometry(data); kw...)
function _zonal(f, x::RasterStackOrArray, ::Nothing, data::Dimension{<:GeometryLookup};
    progress=true, threaded=true, geometrycolumn=nothing, spatialslices, kw...
)
    geoms = data.val.data
    # TODO: filter geoms by raster extent + tree first so we don't descend the
    # full pipeline for geometries outside the raster.

    # If the lookup carries non-default spatial dims and `x` has matching dims
    # by base type, use those instead of the default X/Y.
    if istrue(spatialslices) && data.val.dims != (X(), Y()) &&
        dims(x, data.val.dims) != dims(x, (Val{DD.XDim}(), Val{DD.YDim}()))
        spatialslices = dims(x, data.val.dims)
    end

    n = length(geoms)
    n == 0 && return []
    zs, start_index = _alloc_zonal(f, x, geoms, n; spatialslices, kw...)
    if start_index != n + 1
        _run(start_index:n, threaded, progress, "Applying $f to each geometry...") do i
            zs[i] = _zonal(f, x, geoms[i]; spatialslices, kw...)
        end
    end

    return_lookup_dims = if spatialslices isa DD.AllDims
        dims(data, spatialslices)
    elseif istrue(spatialslices)
        dims(data, (Val{DD.XDim}(), Val{DD.YDim}()))
    else # fallback
        (X(), Y())
    end
    # `X()` here means `X(:)`. Rebuild the dims with `:` so we get a neutral
    # materialised dimension.
    return_lookup = rebuild(lookup(data); dims=rebuild.(return_lookup_dims, (:,)))
    return_dimension = rebuild(data, return_lookup)

    if zs isa AbstractVector{<:Union{<:AbstractDimArray,Missing}}
        return _cat_and_rebuild_parent(x, zs, return_dimension)
    elseif zs isa AbstractVector{<:Union{<:AbstractDimStack,Missing}}
        dimarrays = NamedTuple{names(st)}(
            ntuple(length(names(st))) do i
                _cat_and_rebuild_parent(layers(st)[i], (layers(z)[i] for z in zs), return_dimension)
            end
        )
        return rebuild(x; data=dimarrays, dims=(dims(first(zs))..., return_dimension))
    else
        return Raster(zs, (return_dimension,))
    end
end

# Otherwise of is a geom, table or vector
_zonal(f, x::RasterStackOrArray, of; kw...) = _zonal(f, x, GI.trait(of), of; kw...)

_zonal(f, x, ::GI.AbstractFeatureCollectionTrait, fc; kw...) =
    _zonal(f, x, nothing, fc; kw...)
_zonal(f, x::RasterStackOrArray, ::GI.AbstractFeatureTrait, feature; kw...) =
    _zonal(f, x, GI.geometry(feature); kw...)
function _zonal(f, x::AbstractRaster, ::GI.AbstractGeometryTrait, geom;
    skipmissing, emptyval, spatialslices, missingval, kw...
)
    cropped = crop(x; to=geom, touches=true)
    masked = if length(cropped) == 0
        if istrue(skipmissing) && isfalse(spatialslices)
            return missingval
        end
        cropped # don't mask if we know there is nothing to mask, otherwise it errors
    else
        mask(cropped; with=geom, kw...)
    end
    return _maybe_skipmissing_call(_maybe_spatialsliceify(f, spatialslices, missingval), masked, skipmissing, emptyval)
end
function _zonal(f, st::AbstractRasterStack, ::GI.AbstractGeometryTrait, geom;
    skipmissing, emptyval, spatialslices, missingval, kw...
)
    cropped = crop(st; to=geom, touches=true)
    masked = if length(cropped) == 0
        if istrue(skipmissing) && isfalse(spatialslices)
            return map(_ -> missingval, st)
        end
        cropped # don't mask if we know there is nothing to mask, otherwise it errors
    else
        mask(cropped; with=geom, kw...)
    end
    return maplayers(masked) do A
        if length(A) == 0 && (istrue(skipmissing) && isfalse(spatialslices))
            return missingval
        end
        _maybe_skipmissing_call(_maybe_spatialsliceify(f, spatialslices, missingval), A, skipmissing, emptyval)
    end
end
function _zonal(f, x::RasterStackOrArray, ::Nothing, data; 
    progress=true, threaded=true, geometrycolumn=nothing, missingval, kw...
)
    geoms = _get_geometries(data, geometrycolumn)
    n = length(geoms)
    n == 0 && return []
    zs, start_index = _alloc_zonal(f, x, geoms, n; missingval, kw...)
    start_index == n + 1 && return zs
    _run(start_index:n, threaded, progress, "Applying $f to each geometry...") do i
        zs[i] = _zonal(f, x, geoms[i]; missingval, kw...)
    end

    return_dimension = Dim{:Geometry}(axes(zs, 1))
    
    if zs isa AbstractVector{<: Union{<: AbstractDimArray, Missing}}
        return _cat_and_rebuild_parent(x, zs, return_dimension)
    elseif zs isa AbstractVector{<: Union{<: AbstractDimStack, Missing}}
        dimarrays = NamedTuple{names(st)}(
            ntuple(length(names(st))) do i
                _cat_and_rebuild_parent(layers(st)[i], (layers(z)[i] for z in zs), return_dimension)
            end
        )
        return rebuild(x; data = dimarrays, dims = (dims(first(zs))..., return_dimension))
    end
    return zs
end

function _alloc_zonal(f, x, geoms, n; spatialslices = _True(), missingval, kw...)
    n_missing::Int = 0
    if isfalse(spatialslices)
        # Find first non-missing entry and count number of missing entries
        z1 = _zonal(f, x, first(geoms); spatialslices, missingval, kw...)
        for geom in geoms
            z1 = _zonal(f, x, geom; spatialslices, missingval, kw...)
            if !(ismissing(z1) || z1 === missingval)
                break
            end
            n_missing += 1
        end
        zs = Vector{Union{typeof(missingval), typeof(z1)}}(undef, n)
        zs[1:n_missing] .= missingval
        # Exit early when all elements are missing
        if n_missing == n
            return zs, n_missing + 1
        end
        zs[n_missing + 1] = z1
        return zs, n_missing + 1
    else # spatialslices is true, we know we need an output raster
        z1 = _zonal(f, x, first(geoms); spatialslices, missingval, kw...)
        _missing_array = z1
        for geom in geoms
            z1 = _zonal(f, x, geom; spatialslices, missingval, kw...)
            if !all(ismissing, z1) # here, z1 is a raster, so it's fine - but maybe we should have a better check...
                break
            end
            n_missing += 1
        end
        # TODO: just bite the bullet and use map.  alloc_zonal is a mess.
        zs = Vector{Union{typeof(z1), typeof(_missing_array)}}(undef, n)
        for i in 1:n
            zs[i] = _missing_array
        end
        # Exit early when all elements are missing
        if n_missing == n
            return zs, n_missing + 1
        end
        zs[n_missing + 1] = z1
        return zs, n_missing + 1
    end
end

# No emptyval, just skipmissing or not
_maybe_skipmissing_call(f, A, sm, emptyval::NoKW) = istrue(sm) ? f(skipmissing(A)) : f(A)
# Allow for emptyval if the skipmissing iterator is empty
function _maybe_skipmissing_call(f, A, sm, emptyval)
    if istrue(sm)
        itr = skipmissing(A)
        isempty(itr) && return emptyval
        f(skipmissing(A))
    else
        f(A)
    end
end

# Spatial-slice helpers (`_mapspatialslices`, `_SpatialSliceify`,
# `__do_cat_with_last_dim*`, `_cat_and_rebuild_parent`) live in
# `methods/spatial_slice.jl` and are shared with `extract`.
