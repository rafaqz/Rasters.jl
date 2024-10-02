function _great_circle_bearing(lon1::AbstractFloat, lat1::AbstractFloat, lon2::AbstractFloat, lat2::AbstractFloat)
    dLong = lon1 - lon2

    s = cosd(lat2)*sind(dLong)
    c = cosd(lat1)*sind(lat2) - sind(lat1)*cosd(lat2)*cosd(dLong)

    return atan(s, c)
end

## Get the area of a LinearRing with coordinates in radians
# Using Gidard's theorem
function _linearring_area(ring; radius)
    n = GI.npoint(ring)
    area = -(n-3)*pi

    prevpoint = GI.getpoint(ring, n-1)
    point = GI.getpoint(ring, 1)
    
    for i in 2:n
        nextpoint = GI.getpoint(ring, i)

        beta1 = _great_circle_bearing(GI.x(point), GI.y(point), GI.x(prevpoint), GI.y(prevpoint)) 
        beta2 = _great_circle_bearing(GI.x(point), GI.y(point), GI.x(nextpoint), GI.y(nextpoint)) 
        angle = acos(cos(-beta1)*cos(-beta2) + sin(-beta1)*sin(-beta2))
        area += angle

        prevpoint = point
        point = nextpoint
    end
    
    return area*radius^2
end

_area_from_coords(transform, geom; radius) = _area_from_coords(transform, GI.trait(geom), geom; radius)
function _area_from_coords(transform::AG.CoordTransform, ::GI.LinearRingTrait, ring; radius)
    points = map(GI.getpoint(ring)) do p 
        t = AG.transform!(AG.createpoint(p...), transform)
        (GI.x(t), GI.y(t))
    end
    return _linearring_area(GI.LinearRing(points); radius)
end

# For lat-lon projections. Get the area of each latitudinal band, then multiply by the width
function _area_from_lonlat(lon::XDim, lat::YDim; radius)
    two_pi_R2 = 2 * pi * radius * radius
    band_area = broadcast(DD.intervalbounds(lat)) do yb
        two_pi_R2 * (sind(yb[2]) - sind(yb[1]))
    end
    
    broadcast(DD.intervalbounds(lon), band_area') do xb, ba
        abs(xb[2] - xb[1]) / 360 * ba
    end
end

function _spherical_cellarea(dims::Tuple{<:XDim, <:YDim}; radius = 6371008.8)
    # check the dimensions 
    isnothing(crs(dims)) && _no_crs_error()

    areas = if _isdegrees(crs(dims)) # check if need to reproject
        _area_from_lonlat(dims...; radius)
    elseif !isnothing(mappedcrs(dims)) && _isdegrees(mappedcrs(dims))
        _area_from_lonlat(reproject(dims; crs = mappedcrs(dims))...; radius)
    else
        xbnds, ybnds = DD.intervalbounds(dims)
        AG.crs2transform(crs(dims), EPSG(4326), order = :trad) do transform
            [_area_from_coords(
                transform,         
                GI.LinearRing([
                    (xb[1], yb[1]), 
                    (xb[2], yb[1]), 
                    (xb[2], yb[2]), 
                    (xb[1], yb[2]),
                    (xb[1], yb[1])
                ]); 
                radius
                )
                for xb in xbnds, yb in ybnds]
        end
    end
end

# TODO: put these in ArchGDAL
_isgeographic(crs) = _isgeographic(AG.importCRS(crs))
_isgeographic(crs::AG.ISpatialRef) = AG.GDAL.osrisgeographic(crs) |> Bool

_isdegrees(crs) = _isdegrees(AG.importCRS(crs))
function _isdegrees(crs::AG.ISpatialRef)
    _isgeographic(crs) || return false
    pointer = Ref{Cstring}()
    result = AG.GDAL.osrgetangularunits(crs, pointer)
    return unsafe_string(pointer[]) == "degree"
end