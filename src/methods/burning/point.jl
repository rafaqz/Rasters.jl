# TODO use rasterize_points for this instead

# _fill_point!
# Fill a raster with `fill` where points are inside raster pixels
@noinline _fill_point!(x::RasterStackOrArray, geom; kw...) = _fill_point!(x, GI.geomtrait(geom), geom; kw...)
@noinline function _fill_point!(x::RasterStackOrArray, ::GI.GeometryCollectionTrait, geom; kw...)
    _without_mapped_crs(x) do x1
        for geom in GI.getgeom(geom)
            _fill_point!(x, geom; kw...)
        end
    end
    return true
end
@noinline function _fill_point!(x::RasterStackOrArray, ::GI.AbstractGeometryTrait, geom; kw...)
    # Just find which pixels contain the points, and set them to true
    _without_mapped_crs(x) do x1
        for point in GI.getpoint(geom)
            _fill_point!(x, point; kw...)
        end
    end
    return true
end
@noinline function _fill_point!(x::RasterStackOrArray, ::GI.AbstractPointTrait, point;
    fill, atol=nothing, lock=nothing, kw...
)
    dims1 = commondims(x, DEFAULT_POINT_ORDER)
    selectors = map(dims1) do d
        _at_or_contains(d, _dimcoord(d, point), atol)
    end
    # TODO make a check in dimensionaldata that returns the index if it is inbounds
    if hasselection(x, selectors)
        I = dims2indices(dims1, selectors)
        if isnothing(lock)  
            _fill_index!(x, fill, I)
        else
            sector = CartesianIndices(map(i -> i:i, I))
            Base.lock(lock)
            _fill_index!(x, fill, I)
            Base.unlock(lock)
        end
        return true
    else
        return false
    end
end

function _at_or_contains(d, v, atol)
    selector = sampling(d) isa Intervals ? Contains(v) : At(v; atol=atol)
    DD.basetypeof(d)(selector)
end

# Fill Int indices directly
_fill_index!(st::AbstractRasterStack, fill::NamedTuple, I::NTuple{<:Any,Int}) = st[I...] = fill
_fill_index!(A::AbstractRaster, fill, I::NTuple{<:Any,Int}) = A[I...] = fill
_fill_index!(A::AbstractRaster, fill::Function, I::NTuple{<:Any,Int}) = A[I...] = fill(A[I...])

_fill_index!(st::AbstractRasterStack, fill::NamedTuple, I) = 
    map((A, f) -> A[I...] .= Ref(f), st, fill)
_fill_index!(A::AbstractRaster, fill, I) = A[I...] .= Ref(fill)
_fill_index!(A::AbstractRaster, fill::Function, I) = A[I...] .= fill.(view(A, I...))

# Get the GeoInterface coord from a point for a specific Dimension
_dimcoord(::XDim, point) = GI.x(point)
_dimcoord(::YDim, point) = GI.y(point)
_dimcoord(::ZDim, point) = GI.z(point)
