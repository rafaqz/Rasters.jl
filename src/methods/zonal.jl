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
function zonal(f, x::RasterStack; of, skipmissing=true, spatialslices=_False(), missingval=isnothing(missingval(x)) ? missing : missingval(x), kw...)
    # TODO: open currently doesn't work so well for large rasterstacks,
    # we need to fix that before we can go back to this being a single method
    # on `RasterStackOrArray`.
    _zonal(f, _prepare_for_burning(x), of; skipmissing, spatialslices, missingval, kw...)
end
function zonal(f, x::Raster; of, skipmissing=true, spatialslices=_False(), missingval=isnothing(missingval(x)) ? missing : missingval(x), kw...)
    open(x) do xo
        _zonal(f, _prepare_for_burning(xo), of; skipmissing, spatialslices, missingval, kw...)
    end
end

_zonal(f, x::RasterStackOrArray, of::RasterStackOrArray; kw...) = 
    _zonal(f, x, Extents.extent(of); kw...)
_zonal(f, x::RasterStackOrArray, of::DimTuple; kw...) = 
    _zonal(f, x, Extents.extent(of); kw...)
# We don't need to `mask` with an extent, it's square so `crop` will do enough.
_zonal(f, x::Raster, of::Extents.Extent; skipmissing, spatialslices, missingval) = _maybe_skipmissing_call(_maybe_spatialsliceify(f, spatialslices), crop(x; to=of, touches=true), skipmissing)
function _zonal(f, x::RasterStack, ext::Extents.Extent; skipmissing, spatialslices, missingval)
    cropped = crop(x; to=ext, touches=true)
    if length(cropped) == 0 && skipmissing == true
        return map(_ -> missingval, x)
    end
    return maplayers(cropped) do A
        _maybe_skipmissing_call(_maybe_spatialsliceify(f, spatialslices, missingval), A, skipmissing)
    end
end
# Otherwise of is a geom, table or vector
_zonal(f, x::RasterStackOrArray, of; kw...) = _zonal(f, x, GI.trait(of), of; kw...)

_zonal(f, x, ::GI.AbstractFeatureCollectionTrait, fc; kw...) =
    _zonal(f, x, nothing, fc; kw...)
_zonal(f, x::RasterStackOrArray, ::GI.AbstractFeatureTrait, feature; kw...) =
    _zonal(f, x, GI.geometry(feature); kw...)
function _zonal(f, x::AbstractRaster, ::GI.AbstractGeometryTrait, geom; 
    skipmissing, spatialslices, missingval, kw...
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
    return _maybe_skipmissing_call(_maybe_spatialsliceify(f, spatialslices, missingval), masked, skipmissing)
end
function _zonal(f, st::AbstractRasterStack, ::GI.AbstractGeometryTrait, geom; 
    skipmissing, spatialslices, missingval, kw...
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
        _maybe_skipmissing_call(_maybe_spatialsliceify(f, spatialslices, missingval), A, skipmissing)
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
        _cat_and_rebuild_parent(x, zs, return_dimension)
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

# Optionally wrap the input argument in `skipmissing(A)` is `sm` is true.
_maybe_skipmissing_call(f, A, sm) = istrue(sm) ? f(skipmissing(A)) : f(A)

# the only reason we have AbstractDimArray here is to make sure that DD.otherdims is available.
# We could probably get away with just AbstractArray here otherwise.
# The reason this is not just mapslices is because this drops the sliced dimensions automatically, 
# which is what we want.
function _mapspatialslices(f, x::AbstractDimArray; spatialdims = (Val{DD.XDim}(), Val{DD.YDim}()), missingval = missingval(x))
    dimswewant = DD.otherdims(x, spatialdims)
    if isempty(dimswewant)
        return f(x)
    end
    slicedims = rebuild.(dims(x, dimswewant), axes.((x,), dimswewant))
    if any(isempty, DD.dims(x, spatialdims))
        # If any of the spatial dims are empty, we can just return a constant missing array
        # this way we don't construct the dimslices at all...
        missing_array = FillArrays.Fill{Union{typeof(missingval), eltype(x)}, length(dimswewant)}(missingval, length.(dimswewant))
        return rebuild(x; data = missing_array, dims = dimswewant, refdims = ())
    end
    iterator = (rebuild(x; data = d, dims = dims(d)) for d in DD.DimSlices(x; dims = slicedims, drop = true))
    return rebuild(x; data = f.(iterator), dims = dimswewant, refdims = ())
end
# SkipMissingVal and SkipMissing both store the initial value in the `x` property,
# so we can use the same thing to extract it.
function _mapspatialslices(f, s::Union{SkipMissingVal, Base.SkipMissing}; spatialdims = (Val{DD.XDim}(), Val{DD.YDim}()), missingval = missingval(s.x))
    return _mapspatialslices(f ∘ skipmissing, s.x; spatialdims, missingval)
end
    

_maybe_spatialsliceify(f, spatialslices, missingval = missing) = istrue(spatialslices) ? _SpatialSliceify(f, (Val{DD.XDim}(), Val{DD.YDim}()), missingval) : f
_maybe_spatialsliceify(f, spatialslices::DD.AllDims, missingval = missing) = _SpatialSliceify(f, spatialslices, missingval)

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
struct _SpatialSliceify{F, D, M}
    f::F
    dims::D
    missingval::M
end

(r::_SpatialSliceify{F, D, M})(x) where {F, D, M} = _mapspatialslices(r.f, x; spatialdims = r.dims, missingval = r.missingval)

# This is a helper function that concatenates an array of arrays along their last dimension.
# and returns a ConcatDiskArray so that it doesn't allocate at all.\
# Users can always rechunk later.  But this saves us a lot of time when doing datacube ops.
# And the chunk pattern is available in the concat diskarray.
function __do_cat_with_last_dim(input_arrays)
    # This assumes that the input array is a vector of arrays.
    As = Missings.disallowmissing(collect(input_arrays))
    dims = ndims(first(As)) + 1
    sz = ntuple(dims) do i
        i == dims ? length(As) : 1
    end
    cdas = reshape(As, sz)
    backing_array = DiskArrays.ConcatDiskArray(cdas)
   return backing_array
end

function __do_cat_with_last_dim_multidim_version(As)
    # This CANNOT assume that the input array is a vector of arrays.
    new_n_dims = ndims(As) + ndims(first(As))
    sz = ntuple(new_n_dims) do i
        i <= ndims(first(As)) ? 1 : size(As, i-1)
    end
    cdas = reshape(As, sz)
    backing_array = DiskArrays.ConcatDiskArray(cdas)
   return backing_array
end
# This is a wrapper around the helper function that performs the final cat and rebuild, but on 
# a dimarray.
function _cat_and_rebuild_parent(parent, children, newdim)
    backing_array = __do_cat_with_last_dim(children) # see zonal.jl for implementation
    children_dims = dims(first(children))
    final_dims = DD.format((children_dims..., newdim), backing_array)
    return rebuild(parent; data = backing_array, dims = final_dims)
end

precompile(__do_cat_with_last_dim, (Vector{Raster{<: Any, 1}},))
precompile(__do_cat_with_last_dim, (Vector{Raster{<: Any, 2}},))
precompile(__do_cat_with_last_dim, (Vector{Raster{<: Any, 3}},))
precompile(__do_cat_with_last_dim, (Vector{Raster{<: Any, 4}},))