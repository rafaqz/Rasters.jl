
function _zonal(f, x::RasterStackOrArray, ::Nothing, data::GeometryLookup; kw...)
    return _zonal(f, x, nothing, Dim{:Geometry}(data); kw...)
end
function _zonal(f, x::RasterStackOrArray, ::Nothing, data::Dim{Name, <: GeometryLookup}; 
    progress=true, threaded=true, geometrycolumn=nothing, kw...
) where Name
    geoms = data.val.data
    # TODO: deliberately filter geoms based on extent and tree
    # so that we don't waste time calling `mask` on lots of geometries
    # that are outside the extent of the raster
    # but that is for later.
    n = length(geoms)
    n == 0 && return []
    zs, start_index = _alloc_zonal(f, x, geoms, n; kw...)
    start_index == n + 1 && return zs
    _run(start_index:n, threaded, progress, "Applying $f to each geometry...") do i
        zs[i] = _zonal(f, x, geoms[i]; kw...)
    end
    @show typeof(zs)
    if zs isa AbstractVector{<: Union{<: AbstractDimArray, Missing}}
        println("Got a vector of DimArrays")
        return cat(zs...; dims = data)
    else
        println("Got a vector of other things")
        return Raster(zs, (data,))
    end
end
