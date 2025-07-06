#=
## Reproject

Reproject just forwards to `GO.reproject`.  Since regular reproject here in Rasters and GO reproject both need Proj,
this is not too bad.
=#

function reproject(target::GeoFormat, l::GeometryLookup)
    source = l.crs
    # TODO: allow GDAL reproject for its antimeridian cutting
    return rebuild(l; data = GO.reproject(l.data; source_crs = source, target_crs = target, always_xy = true), crs = target)
end

#=
## Zonal

Zonal with a geometry lookup or a geometry dimension should return a vector data cube.
=#
function _zonal(f, x::RasterStackOrArray, ::Nothing, data::GeometryLookup; kw...)
    return _zonal(f, x, nothing, Geometry(data); kw...)
end
function _zonal(f, x::RasterStackOrArray, ::Nothing, data::Dimension{<: GeometryLookup}; 
    progress=true, threaded=true, geometrycolumn=nothing, spatialslices, kw...
)
    geoms = data.val.data
    # TODO: deliberately filter geoms based on extent and tree
    # so that we don't waste time descending through the pipeline
    # for geometries that are outside the extent of the raster
    # but that is for later.

    if istrue(spatialslices) && data.val.dims != (X(), Y()) && dims(x, data.val.dims) != dims(x, (Val{DD.XDim}(), Val{DD.YDim}()))
        spatialslices = dims(x, data.val.dims) # use the dimensions in the geometry lookup!
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
    # Note here that e.g. `X()` is actually `X(:)`.  That's why we rebuild the dims with colons - 
    # so that we get a "neutral materialized dimension" out of it.
    return_lookup = rebuild(lookup(data); dims = rebuild.(return_lookup_dims, (:,)))

    return_dimension = rebuild(data, return_lookup)

    if zs isa AbstractVector{<: Union{<: AbstractDimArray, Missing}}
        return _cat_and_rebuild_parent(x, zs, return_dimension)
    elseif zs isa AbstractVector{<: Union{<: AbstractDimStack, Missing}}
        dimarrays = NamedTuple{names(st)}(
            ntuple(length(names(st))) do i
                _cat_and_rebuild_parent(layers(st)[i], (layers(z)[i] for z in zs), return_dimension)
            end
        )
        return rebuild(x; data = dimarrays, dims = (dims(first(zs))..., return_dimension))
    else
        return Raster(zs, (return_dimension,))
    end
end

#=
## Crop

Cropping a geometry lookup to either a geometry lookup, a geometry, or a bounding box should return a geometry lookup.
=#
