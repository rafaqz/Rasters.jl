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
    # don't assume the ring is a GI ring
    (p1, p2, p3, p4) = _lonlat_to_sphericalpoint.(ring)
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

# _area_from_coords(transform, geom) = _area_from_coords(transform, geom)
function _area_from_coords(transform::Proj.Transformation, ring_points)
    t = (transform(GI.x(p), GI.y(p)) for p in ring_points)
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

function _spherical_cellarea(dims::Tuple{<:XDim, <:YDim}; radius = 6371008.8, use_area_of_use = true)
    # check the dimensions 
    isnothing(crs(dims)) && _no_crs_error()

    areas = if _isdegrees(crs(dims)) # check if need to reproject
        _area_from_lonlat(dims...; radius)
    elseif !isnothing(mappedcrs(dims)) && _isdegrees(mappedcrs(dims))
        _area_from_lonlat(reproject(dims; crs = mappedcrs(dims))...; radius)
    else
        xbnds, ybnds = DD.intervalbounds(dims)
        R2 = radius * radius
        area_of_use = if use_area_of_use
            # Use a temporary transformation object to get the area of use,
            # since the transformation object has to be recreated with the area of use
            _get_area_of_use(Proj.Transformation(crs(dims), EPSG(4326); always_xy = true), Extents.Extent(X = extrema(dims[1]), Y = extrema(dims[2])))
        else
            C_NULL
        end

        transform = Proj.Transformation(crs(dims), EPSG(4326); always_xy = true, area = area_of_use)

        result = [
                _area_from_coords(
                    transform,         
                    (
                        (xb[1], yb[1]), 
                        (xb[2], yb[1]), 
                        (xb[2], yb[2]), 
                        (xb[1], yb[2]),
                    )
                ) * R2
                for xb in xbnds, yb in ybnds
            ]
        
        if use_area_of_use
            Proj.proj_area_destroy(area_of_use)
        end

        return result
    end
end

function _get_area_of_use(transform::Proj.Transformation, extent::Extents.Extent; densify_pts = 21)
    # Transform the extent using `proj_trans_bounds`
    (xmin, xmax), (ymin, ymax) = Proj.bounds(transform, extent.X, extent.Y; densify_pts)

    # Create an area of use object 
    # This MUST be destroyed by the caller
    area = Proj.proj_area_create()

    # Set the bounding box of the area of use
    Proj.proj_area_set_bbox(
        area, 
        xmin, # west_lon_degree
        ymin, # south_lat_degree
        xmax, # east_lon_degree
        ymax  # north_lat_degree
    )

    return area
end

# TODO: put these in Proj (specifically the dispatches on GFT types)
_isgeographic(crs) = _isgeographic(Proj.CRS(crs))
_isgeographic(crs::Proj.CRS) = Proj.is_geographic(crs)

_isdegrees(crs) = _isdegrees(Proj.CRS(crs))
function _isdegrees(crs::Proj.CRS)
    _isgeographic(crs) || return false
    # This is a tiny bit inefficient, but it takes 500ns on my machine, 
    # so I think we can disregard the inefficiency...
    return axis_is_degrees(crs, 0) && axis_is_degrees(crs, 1)
end

function axis_is_degrees(crs::Proj.CRS, axis_index::Int; context::Ptr{Proj.PJ_CONTEXT} = C_NULL)
    @assert axis_index in (0, 1)
    # In Proj, the `CoordinateSystem` object is contained within the `CRS` object.
    # So we need to extract it explicitly, since Proj doesn't provide utilities for this.
    cs = Proj.proj_crs_get_coordinate_system(crs.pj, context)

    # Instantiate some refs that we'll use to get information out of the PJ struct
    auth_name = Ref{Cstring}()
    code = Ref{Cstring}()
    unitname = Ref{Cstring}()
    
    # Load unit info for the given axis into the pointers
    Proj.proj_cs_get_axis_info(
        cs, 
        axis_index, 
        C_NULL, # out_name
        C_NULL, # out_abbrev
        C_NULL, # out_direction
        C_NULL, # out_unit_conv_factor
        C_NULL, # out_unit_name
        auth_name,
        code,
        context
    )

    # We don't `unsafe_string` the C strings, because we're just going to pass them to Proj's unit lookup function
    Proj.proj_uom_get_info_from_database(
        auth_name[],
        code[],
        unitname,   # out_name
        C_NULL,     # out_conv_factor
        C_NULL,     # out_category
        context
    )
    # Destroy the coordinate system object
    Proj.proj_destroy(cs)

    unit_str = unsafe_string(unitname[])

    # TODO: check if this is the correct way to check if the unit is degrees
    # We can also check if the unit category is "angular", but I chose to replicate
    # the original test from ArchGDAL here.
    # If the unit is not "degree", we could still have an angular unit (radians or some linearly scaled thing),
    # in which case we should technically return true.
    # We'd also have to return the conversion factor in this case, and maybe the category (radians or degrees)...
    return isequal(unit_str, "degree")
end