## Get the area of a LinearRing with coordinates in radians
struct SphericalPoint{T <: Real}
	data::NTuple{3, T}
end
SphericalPoint(x, y, z) = SphericalPoint((x, y, z))

# define the 4 basic mathematical operators elementwise on the data tuple
Base.:+(p::SphericalPoint, q::SphericalPoint) = SphericalPoint(p.data .+ q.data)
Base.:-(p::SphericalPoint, q::SphericalPoint) = SphericalPoint(p.data .- q.data)
Base.:*(p::SphericalPoint, q::SphericalPoint) = SphericalPoint(p.data .* q.data)
Base.:/(p::SphericalPoint, q::SphericalPoint) = SphericalPoint(p.data ./ q.data)
# Define sum on a SphericalPoint to sum across its data
Base.sum(p::SphericalPoint) = sum(p.data)

# define dot and cross products
dot(p::SphericalPoint, q::SphericalPoint) = sum(p * q)
function cross(a::SphericalPoint, b::SphericalPoint)
	a1, a2, a3 = a.data
    b1, b2, b3 = b.data
	SphericalPoint((a2*b3-a3*b2, a3*b1-a1*b3, a1*b2-a2*b1))
end

function _spherical_quadrilateral_area(ring)
    (p1, p2, p3, p4) = _lonlat_to_sphericalpoint.(GI.getpoint(ring))
    area = 0.0
    area += _spherical_triangle_area(p1, p2, p3)
    area += _spherical_triangle_area(p3, p4, p1)
end

# Using Eriksson's formula for the area of spherical triangles: https://www.jstor.org/stable/2691141
function _spherical_triangle_area(a, b, c)
    #t = abs(dot(a, cross(b, c)))
    #t /= 1 + dot(b,c) + dot(c, a) + dot(a, b)
    t = abs(dot(a, (cross(b - a, c - a))) / dot(b + a, c + a))
    2*atan(t)
end

_lonlat_to_sphericalpoint(args) = _lonlat_to_sphericalpoint(args...)
function _lonlat_to_sphericalpoint(lon, lat)
    x = cosd(lat) * cosd(lon)
    y = cosd(lat) * sind(lon)
    z = sind(lat)
    return SphericalPoint(x,y,z)
end

_area_from_coords(transform, geom) = _area_from_coords(transform, GI.trait(geom), geom)
function _area_from_coords(transform::AG.CoordTransform, trait::GI.AbstractCurveTrait, ring)
    t = AG.transform!(GI.convert(AG.geointerface_geomtype(trait), ring), transform)
    return _spherical_quadrilateral_area(t)
end

# For lat-lon projections. Get the area of each latitudinal band, then multiply by the width
function _area_from_lonlat(lon::XDim, lat::YDim; radius)
    two_pi_R2 = 2 * pi * radius * radius
    band_area = broadcast(DD.intervalbounds(lat)) do yb
        two_pi_R2 * (sind(yb[2]) - sind(yb[1]))
    end
    
    broadcast(DD.intervalbounds(lon), band_area') do xb, ba
        abs((xb[2] - xb[1]) / 360 * ba)
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
        R2 = radius * radius
        AG.crs2transform(crs(dims), EPSG(4326); order = :trad) do transform
            [_area_from_coords(
                transform,         
                GI.LinearRing([
                    (xb[1], yb[1]), 
                    (xb[2], yb[1]), 
                    (xb[2], yb[2]), 
                    (xb[1], yb[2]),
                ])
                ) * R2
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