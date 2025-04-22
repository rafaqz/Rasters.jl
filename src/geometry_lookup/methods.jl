
function _zonal(f, x::RasterStackOrArray, ::Nothing, data::GeometryLookup; kw...)
    return _zonal(f, x, nothing, Geometry(data); kw...)
end
function _zonal(f, x::RasterStackOrArray, ::Nothing, data::Dimension{<: GeometryLookup}; 
    progress=true, threaded=true, geometrycolumn=nothing, spatialslices, kw...
)
    geoms = data.val.data
    # TODO: deliberately filter geoms based on extent and tree
    # so that we don't waste time calling `mask` on lots of geometries
    # that are outside the extent of the raster
    # but that is for later.

    if istrue(spatialslices) && data.val.dims != (X(), Y()) && dims(x, data.val.dims) != dims(x, (Val{DD.XDim}(), Val{DD.YDim}()))
        spatialslices = dims(x, data.val.dims) # use the dimensions in the geometry lookup!
    end

    n = length(geoms)
    n == 0 && return []
    zs, start_index = _alloc_zonal(f, x, geoms, n; spatialslices, kw...)
    start_index == n + 1 && return zs
    _run(start_index:n, threaded, progress, "Applying $f to each geometry...") do i
        zs[i] = _zonal(f, x, geoms[i]; spatialslices, kw...)
    end

    return_lookup_dims = if istrue(spatialslices)
        dims(data, (Val{DD.XDim}(), Val{DD.YDim}()))
    elseif spatialslices isa DD.AllDims
        dims(data, spatialslices)
    else # fallback
        (X(), Y())
    end

    return_lookup = rebuild(data.val; dims = rebuild.(return_lookup_dims, (:,)))

    return_dimension = rebuild(data, return_lookup)

    if zs isa AbstractVector{<: Union{<: AbstractDimArray, <: AbstractDimStack, Missing}}
        backing_array = __do_cat_with_last_dim(zs)
        z_dims = dims(first(zs))
        new_dims = DD.format((z_dims..., return_dimension), backing_array)
        return rebuild(x; data = backing_array, dims = new_dims)
    else
        # TODO: how should we reconstruct a rasterstack from a vector of named tuples?
        return Raster(zs, (return_dimension,))
    end
end