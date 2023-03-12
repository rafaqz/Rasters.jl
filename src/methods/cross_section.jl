
"""
    cross_section(x::RasterStackOrArray; with)

Returns the cross section of a `Raster` or `RasterStack` and a line.

The new dimension is call `:cross_section`.

All geometries can be used, but they will all be treated as lines.

# Arguments

- `x`: A `Raster` or `RasterStack`
- `with`: A GeoInterface.jl compatible geometry, Tables.jl compatible table with a `:geometry` 
    column, or an `AbstractVector` of any of these objects..

# Example

```julia
using Rasters, Shapefile, Downloads, Statistics, DataFrames

# Download a borders shapefile
shp_url = "https://github.com/nvkelso/natural-earth-vector/raw/master/10m_cultural/ne_10m_admin_0_countries.shp"
dbf_url = "https://github.com/nvkelso/natural-earth-vector/raw/master/10m_cultural/ne_10m_admin_0_countries.dbf"
shp_name, dbf_name = "country_borders.shp", "country_borders.dbf"
isfile(shp_name) || Downloads.download(shp_url, shp_name)
isfile(dbf_url) || Downloads.download(dbf_url, dbf_name)

# Download and read a raster stack
st = RasterStack(WorldClim{Climate}; month=Jan, lazy=false)


```
"""
cross_section(x::RasterStackOrArray; with) = _cross_section(x, GI.trait(with), with)

_cross_section(x::RasterStackOrArray, with) = _cross_section(x, GI.trait(with), with)
function _cross_section(x, ::GI.AbstractFeatureCollectionTrait, fc)
    [cross_section(x, feature) for feature in GI.getfeature(fc)]
end
_cross_section(x, ::GI.AbstractFeatureTrait, obj) = _cross_section(x, GI.geometry(obj))
function _cross_section(x, ::Nothing, obj)
    if Tables.istable(obj)
        geoms = Tables.getcolumn(obj, first(GI.geometrycolumns(obj)))
        return [cross_section(x, geom) for geom in geoms]
    elseif obj isa AbstractVector
        return [cross_section(x, geom) for geom in obj]
    else
        throw(ArgumentError("Cannot calculate cross_section statistics for objects of type $(typeof(obj))"))
    end
end
function _cross_section(obj::RasterStackOrArray, ::GI.AbstractGeometryTrait, geom)
    xd, yd = xydims = dims(obj, DEFAULT_POINT_ORDER)
    inds = Tuple{Int,Int}[]
    # Combine vectors of indices of each line segment
    _with_lines(geom) do i, line
        append!(inds, _line_inds(xydims, line))
    end
    cross_dim = _cross_dim(obj, inds, xydims)
    # Take a view for each index
    if length(dims(obj)) > 2
        # Concat a series of arrays/stacks
        _cross_section_views(obj, inds, xydims, cross_dim)
    else
        # get values
        if obj isa AbstractRasterStack
            map(NamedTuple(obj)) do A
                _cross_section_values(A, inds, xydims, cross_dim)
            end
        else
            _cross_section_values(obj, inds, xydims, cross_dim)
        end
    end
end

function _cross_section_values(A, inds, xydims, cross_dim)
    cells = map(inds) do I
        D = map((d, i) -> basetypeof(d)(i); xydims, I)
        A[D...]
    end
    Raster(cells, (cross_dim,))
end

function _cross_section_views(A, inds, xydims, cross_dim)
    cells = map(inds) do I
        D = map(xydims, I) do d, i
            basetypeof(d)(i)
        end
        Base.view(A, D...)
    end
    return combine(RasterSeries(cells, cross_dim), cross_dim)
end

function _cross_dim(x, inds, xydims)
    all_points = DimPoints(dims(x, xydims))
    cs_points = map(inds) do I
        D = map(xydims, I) do d, i
            basetypeof(d)(i)
        end
        all_points[dims2indices(dims(all_points), D)...]
    end
    return Dim{:cross_section}(Dimensions.MergeLookup(cs_points, map(DD.basedims, xydims)))
end

