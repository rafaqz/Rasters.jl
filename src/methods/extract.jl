"""
   extract(x, points; order, atol)

Extracts the value of `Raster` or `RasterStack` at given points, returning
a vector of `NamedTuple` with columns for the point geometry and values
value/s.

Note that if objects have more dimensions than the length of the point tuples,
sliced arrays or stacks will be returned instead of single values.

# Arguments

- `x`: a `Raster` or `RasterStack` to extract values from.
- `points`: multiple `Vector`s of point values, a `Vector{Tuple}`,
    or a single `Tuple` or `Vector`. `points` can also be a Tables.jl compatible
    table, in which case `order` may need to specify the keys.

# Keywords

- `order`: a tuple of `Dimension` connecting the order of the points to the array
    axes, such as `(X, Y)`, with a defaut `(XDim, YDim, ZDim)` order.
    If `points` is a table, `order` should be a `Tuple` of `Dimension`/`Symbol` pairs
    like `(X => :xcol, Y => :ycol)`. This will be automatically detected wherever
    possible, assuming the keys match the dimensions of the object `x`.
- `atol`: a tolorerance for floating point lookup values for when the `LookupArray`
    contains `Points`. `atol` is ignored for `Intervals`.

Note: extracting polygons in a `GeoInterface.AbstractGeometry` is not yet supported,
but will be in future.

# Example

Here we extact points matching the occurrence of the Mountain Pygmy Possum,
_Burramis parvus_. This could be used to fit a species distribution lookupl.

```jldoctest
using Rasters, GBIF, CSV

# Get a stack of BioClim layers, and replace missing values with `missing`
st = RasterStack(WorldClim{BioClim}, (1, 3, 5, 7, 12))[Band(1)] |> replace_missing

# Download some occurrence data
obs = GBIF.occurrences("scientificName" => "Burramys parvus", "limit" => 5)

# use `extract` to get values for all layers at each observation point.
points = map(o -> (o.longitude, o.latitude), obs)
vals = extract(st, points)

# output
5-element Vector{NamedTuple{(:X, :Y, :bio1, :bio3, :bio5, :bio7, :bio12)}}:
 (X = missing, Y = missing, bio1 = missing, bio3 = missing, bio5 = missing, bio7 = missing, bio12 = missing)
 (X = 147.096394, Y = -36.935687, bio1 = 9.408354f0, bio3 = 40.790546f0, bio5 = 22.39425f0, bio7 = 23.0895f0, bio12 = 1292.0f0)
 (X = 148.450743, Y = -35.999643, bio1 = 8.269542f0, bio3 = 41.030262f0, bio5 = 21.4485f0, bio7 = 23.858f0, bio12 = 1440.0f0)
 (X = 148.461854, Y = -36.009001, bio1 = 6.928167f0, bio3 = 41.78015f0, bio5 = 20.18025f0, bio7 = 23.69975f0, bio12 = 1647.0f0)
 (X = 148.459452, Y = -36.002648, bio1 = 6.928167f0, bio3 = 41.78015f0, bio5 = 20.18025f0, bio7 = 23.69975f0, bio12 = 1647.0f0)

```
"""
function extract(x::RasterStackOrArray, data; 
    dims=DD.dims(x, DEFAULT_POINT_ORDER), kw...
)
    _extract(x, data; dims, names=_names(x), kw...)
end
_extract(A::RasterStackOrArray, point::Missing; kw...) = missing
function _extract(A::RasterStackOrArray, geom; kw...) 
    _extract(GI.geomtrait(geom), A, geom; kw...) 
end
function _extract(::GI.AbstractFeatureTrait, A::RasterStackOrArray, feature; kw...) 
    _extract(A, GI.geometry(feature); kw...) 
end
function _extract(::GI.AbstractMultiPointTrait, A::RasterStackOrArray, geom; kw...) 
    (_extract(A, p; kw...) for p in GI.getpoint(geom))
end
function _extract(::Nothing, A::RasterStackOrArray, geoms; kw...) 
    GI.isgeometry(first(geoms))
    (_extract(A, g; kw...) for g in geoms)
end
function _extract(::GI.AbstractGeometryTrait, A::RasterStackOrArray, geom; names, kw...) 
    B = boolmask(geom; to=dims(A, DEFAULT_POINT_ORDER), kw...)
    pts = DimPoints(B)
    dis = DimIndices(B)
    ((; geometry=_geom_nt(dims(B), pts[I]), _prop_nt(A, I, names)...) for I in CartesianIndices(B) if B[I])  
end
_geom_nt(dims::DimTuple, pts) = NamedTuple{map(dim2key, dims)}(pts) 
_prop_nt(st::AbstractRasterStack, I, names::NamedTuple{K}) where K = NamedTuple{K}(values(st[I]))
_prop_nt(A::AbstractRaster, I, names::NamedTuple{K}) where K = NamedTuple{K}((A[I],))

function _extract(::GI.PointTrait, x::RasterStackOrArray, point; dims, names, atol=nothing)
    # Get the actual dimensions available in the object
    coords = map(DD.dims(x)) do d
        _dimcoord(d, point)
    end

    # Extract the values
    if any(map(ismissing, coords)) 
        point_vals = map(_ -> missing, coords)
        layer_vals = map(_ -> missing, layer_keys)
    else
        selectors = map(dims, coords) do d, c
            _at_or_contains(d, c, atol)
        end
        layer_vals = if DD.hasselection(x, selectors)
            x isa Raster ? (x[selectors...],) : x[selectors...]
        else
            map(_ -> missing, names)
        end
    end
    geometry = point
    properties = NamedTuple{keys(names)}(layer_vals)
    return (; geometry, properties...)
end

_names(A::AbstractRaster) = NamedTuple{(Symbol(name(A)),)}((Symbol(name(A)),))
_names(A::AbstractRasterStack) = NamedTuple{keys(A)}(keys(A))
