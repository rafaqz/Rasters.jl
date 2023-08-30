"""
    cellsize(x)

Gives the approximate size of each cell in square km. 
This function works for any projection, but uses an algorithm for spherical polygons, and thus gives the approximates the true size to about 0.1%, depending on latitude.

# Arguments

- `x`: A `Raster` or a `Tuple` of `X` and `Y` dimensions. 

## Example

```jldoctest
using Rasters, Rasters.LookupArrays
dimz = X(Projected(90.0:0.1:120; sampling=Intervals(Start()), order = ForwardOrdered(), span = Regular(0.1), crs=EPSG(4326))),
       Y(Projected(0.0:0.1:50; sampling=Intervals(Start()), order = ForwardOrdered(), span = Regular(0.1), crs=EPSG(4326)))

cs = cellsize(dimz)

```
$EXPERIMENTAL
"""

function _great_circle_bearing(lon1::AbstractFloat, lat1::AbstractFloat, lon2::AbstractFloat, lat2::AbstractFloat)
    dLong = lon1 - lon2

    s = cos(lat2)*sin(dLong)
    c = cos(lat1)*sin(lat2) - sin(lat1)*cos(lat2)*cos(dLong)

    return atan(s, c)
end

## Get the area of a LinearRing with coordinates in radians
# Using Gidard's theorem
function _area_from_rads(ring; R = 6371.0088)
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
    
    return area*R^2
end

_area_from_coords(transform, geom) = _area_from_coords(transform, GI.trait(geom), geom)
_area_from_coords(geom) = _area_from_coords(GI.trait(geom), geom)

function _area_from_coords(::GI.LinearRingTrait, ring) # no ArchGDAL, assumes degrees
    points = map(GI.getpoint(ring)) do p 
        (deg2rad(GI.x(p)), deg2rad(GI.y(p)))
    end

    return _area_from_rads(GI.LinearRing(points))
end

function _area_from_coords(transform::ArchGDAL.CoordTransform, ::GI.LinearRingTrait, ring)
    points = map(GI.getpoint(ring)) do p 
        t = ArchGDAL.transform!(ArchGDAL.createpoint(p...), transform)
        (deg2rad(GI.x(t)), deg2rad(GI.y(t)))
    end
    return _area_from_rads(GI.LinearRing(points))
end

function cellsize(x::Tuple{X, Y})
    xbnds, ybnds = DimensionalData.intervalbounds(x)

    if convert(CoordSys, crs(x)) == CoordSys("Earth Projection 1, 104") # check if need to reproject
        areas = [_area_from_coords(
            GI.LinearRing([
                (xb[1], yb[1]), 
                (xb[2], yb[1]), 
                (xb[2], yb[2]), 
                (xb[1], yb[2]),
                (xb[1], yb[1])
            ]))
            for xb in xbnds, yb in ybnds]
    else 
        areas = ArchGDAL.crs2transform(crs(x), EPSG(4326)) do transform
            [_area_from_coords(
                transform,         
                GI.LinearRing([
                    (xb[1], yb[1]), 
                    (xb[2], yb[1]), 
                    (xb[2], yb[2]), 
                    (xb[1], yb[2]),
                    (xb[1], yb[1])
                ]))
                for xb in xbnds, yb in ybnds]
        end
    end

    return Raster(areas, dims = x)
end

function cellsize(x::Raster)
    cellsize(dims(x, (X, Y)))
end
