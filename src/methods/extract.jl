"""
   extract(x, points; order, atol)

Extracts the value of `Raster` or `RasterStack` at given points, returning
a vector of `NamedTuple` with columns for the point dimensions and layer
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
function extract(A::RasterStackOrArray, points::NTuple{<:Any,<:AbstractVector}; kw...)
    extract(A, zip(points...); kw...)
end
function extract(A::RasterStackOrArray, points::GI.AbstractGeometry; kw...)
    extract(A, _flat_nodes(GI.coordinates(points)); kw...)
end
function extract(A::RasterStackOrArray, points::AbstractVector{<:Tuple}; kw...)
    extract.(Ref(A), points; kw...)
end
function extract(A::RasterStackOrArray, points::AbstractVector{Union{<:AbstractVector{<:Union{Real,Missing}},Missing}}; kw...)
    extract.(Ref(A), points; kw...)
end
function extract(A::RasterStackOrArray, data; order=nothing, kw...) 
    if Tables.istable(data)
        order = isnothing(order) ? _auto_dim_columns(dims(A), data) : order
        rows = Tables.rows(data)
        point_dims = map(p -> DD.basetypeof(p[1])(p[2]), order)
        point_keys = map(val, point_dims)
        map(rows) do row
            point_vals = map(pk -> row[pk], point_keys)
            extract(A, point_vals; order=map(first, order), point_keys, kw...)
        end
    else
        order = isnothing(order) ? DEFAULT_POINT_ORDER : order
        map(data) do point
            extract(A, point; order, kw...)
        end
    end
end
extract(A::RasterStackOrArray, points::Missing; kw...) = missing
function extract(
    x::RasterStackOrArray, point::Union{Tuple,AbstractVector{<:Union{Missing,<:Real}}};
    order=(XDim, YDim, ZDim),
    point_keys=map(DD.dim2key, dims(x, order)),
    layer_keys=_layer_keys(x),
    atol=nothing
)
    # Get the actual dimensions available in the object
    # Usually this will be `X` and `Y`, but `Z` as well if it exists.
    ordered_dims = dims(x, order)
    length(point) == length(ordered_dims) || throw(ArgumentError("Length of `point` does not match dims. Pass `order` dims manually"))
    point = ntuple(i -> point[i], length(ordered_dims))
    dimtypes = map(DD.basetypeof, ordered_dims)

    # Extract the values
    if any(map(ismissing, point)) 
        point_vals = map(_ -> missing, ordered_dims)
        layer_vals = map(_ -> missing, layer_keys)
    else
        selectors = map((d, x) -> _at_or_contains(d, x, atol), ordered_dims, point)
        point_vals = map(val ∘ val, selectors)
        layer_vals = if DD.hasselection(x, selectors)
            x isa Raster ? (x[selectors...],) : x[selectors...]
        else
            map(_ -> missing, layer_keys)
        end
    end
    return NamedTuple{(point_keys..., layer_keys...)}((point_vals..., layer_vals...))
end


_layer_keys(A::AbstractRaster) = cleankeys(name(A))
_layer_keys(A::AbstractRasterStack) = cleankeys(keys(A))
