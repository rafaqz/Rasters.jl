const XYExtent = Extents.Extent{(:X,:Y),Tuple{Tuple{Float64,Float64},Tuple{Float64,Float64}}}

# Get the bounds of a geometry
_extent(geom; kw...)::XYExtent = _extent(GI.trait(geom), geom; kw...)
function _extent(::Nothing, data::AbstractVector; kw...)::XYExtent
    g1 = first(data)
    if GI.trait(g1) isa GI.PointTrait 
        xs = extrema(p -> GI.x(p), data)
        ys = extrema(p -> GI.y(p), data)
        return _float64_xy_extent(Extents.Extent(X=xs, Y=ys))
    else
        ext = reduce(data; init=_extent(first(data))) do ext, geom
            Extents.union(ext, _extent(geom))
        end
        return _float64_xy_extent(ext)
    end
end
_extent(::Nothing, data::RasterStackOrArray; kw...)::XYExtent = _float64_xy_extent(Extents.extent(data))
function _extent(::Nothing, data::T; geometrycolumn=nothing)::XYExtent where T
    if Tables.istable(T)
        singlecolumn = isnothing(geometrycolumn) ? first(GI.geometrycolumns(data)) : geometrycolumn
        cols = Tables.columns(data)
        if singlecolumn isa Symbol && singlecolumn in Tables.columnnames(cols)
            # Table of geometries
            geoms = Tables.getcolumn(data, singlecolumn)
            return _extent(nothing, geoms)
        else
            multicolumn = isnothing(geometrycolumn) ? DEFAULT_POINT_ORDER : geometrycolumn 
            # TODO: test this branch
            # Table of points with dimension columns
            bounds = reduce(multicolumn; init=(;)) do acc, key
                if key in Tables.columnnames(cols)
                    merge(acc, (; key=extrema(cols[key])))
                else
                    acc
                end
            end
            return _float64_xy_extent(Extents.Extent(bounds))
        end
    else
        ext = Extents.extent(data)
        ext isa Extents.Extent || throw(ArgumentError("object returns `nothing` from `Extents.extent`."))
        return _float64_xy_extent(ext)
    end
end
function _extent(::GI.AbstractPointTrait, point; kw...)::XYExtent
    x, y = Float64(GI.x(point)), Float64(GI.y(point))
    Extents.Extent(X=(x, x), Y=(y, y))
end
function _extent(::GI.AbstractGeometryTrait, geom; kw...)::XYExtent
    geomextent = GI.extent(geom; fallback=false)
    if isnothing(geomextent)
        points = GI.getpoint(geom)
        xbounds = extrema(GI.x(p) for p in points)
        ybounds = extrema(GI.y(p) for p in points)
        return _float64_xy_extent(Extents.Extent(X=xbounds, Y=ybounds))
    else
        return _float64_xy_extent(geomextent)
    end
end
_extent(::GI.AbstractFeatureTrait, feature; kw...)::XYExtent = _extent(GI.geometry(feature))
function _extent(::GI.GeometryCollectionTrait, collection; kw...)::XYExtent
    geometries = GI.getgeom(collection)
    init = _float64_xy_extent(_extent(first(geometries)))
    ext = reduce(geometries; init) do acc, g
        Extents.union(acc, _extent(g))
    end
    return _float64_xy_extent(ext)
end
function _extent(::GI.AbstractFeatureCollectionTrait, features; kw...)::XYExtent
    features = GI.getfeature(features)
    init = _float64_xy_extent(_extent(first(features)))
    ext = reduce(features; init) do acc, f
        Extents.union(acc, _extent(f))
    end
    return _float64_xy_extent(ext)
end

function _float64_xy_extent(ext::Extents.Extent)
    xbounds = map(Float64, ext.X)
    ybounds = map(Float64, ext.Y)
    return Extents.Extent(X=xbounds, Y=ybounds)
end
