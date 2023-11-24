"""
    extract(x, geoms; atol)

Extracts the value of `Raster` or `RasterStack` at given points, returning
an iterable of `NamedTuple` with properties for `:geometry` and raster or
stack layer values.

Note that if objects have more dimensions than the length of the point tuples,
sliced arrays or stacks will be returned instead of single values.

# Arguments

- `x`: a `Raster` or `RasterStack` to extract values from.
- `geoms`: GeoInterface.jl compatible geometries, or tables or iterables of geometries.

# Keywords

- `geometry`: include a `:geometry` column with the corresponding points for each value, `true` by default.
- `index`: include a column of the `CartesianIndex` for each value, `false` by default.
- `names`: `Tuple` of `Symbol` corresponding to layers of a `RasterStack`. All layers by default.
- `skipmissing`: skip missing points automatically.
- `atol`: a tolorerance for floating point lookup values for when the `LookupArray`
    contains `Points`. `atol` is ignored for `Intervals`.

# Example

Here we extact points matching the occurrence of the Mountain Pygmy Possum,
_Burramis parvus_. This could be used to fit a species distribution model.

```julia
using Rasters, RasterDataSources, ArchGDAL, GBIF2, CSV

# Get a stack of BioClim layers, and replace missing values with `missing`
st = RasterStack(WorldClim{BioClim}, (1, 3, 5, 7, 12)) |> replace_missing

# Download some occurrence data
obs = GBIF2.occurrence_search("Burramys parvus"; limit=5, year="2009")

# use `extract` to get values for all layers at each observation point.
# We `collect` to get a `Vector` from the lazy iterator.
collect(extract(st, pnts; skipmissing=true))

# output
5-element Vector{NamedTuple{(:geometry, :bio1, :bio3, :bio5, :bio7, :bio12)}}:
 (geometry = (0.21, 40.07), bio1 = 17.077084f0, bio3 = 41.20417f0, bio5 = 30.1f0, bio7 = 24.775f0, bio12 = 446.0f0)
 (geometry = (0.03, 39.97), bio1 = 17.076923f0, bio3 = 39.7983f0, bio5 = 29.638462f0, bio7 = 24.153847f0, bio12 = 441.0f0)
 (geometry = (0.03, 39.97), bio1 = 17.076923f0, bio3 = 39.7983f0, bio5 = 29.638462f0, bio7 = 24.153847f0, bio12 = 441.0f0)
 (geometry = (0.52, 40.37), bio1 = missing, bio3 = missing, bio5 = missing, bio7 = missing, bio12 = missing)
 (geometry = (0.32, 40.24), bio1 = 16.321388f0, bio3 = 41.659454f0, bio5 = 30.029825f0, bio7 = 25.544561f0, bio12 = 480.0f0)
```
"""
function extract end
function extract(x::RasterStackOrArray, data;
    dims=DD.dims(x, DEFAULT_POINT_ORDER), names=_names(x), kw...
)
    _extract(x, data; dims, names, kw...)
end

_extract(A::RasterStackOrArray, point::Missing; kw...) = missing
_extract(A::RasterStackOrArray, geom; kw...) = _extract(A, GI.geomtrait(geom), geom; kw...)
function _extract(A::RasterStackOrArray, ::Nothing, geoms; skipmissing=false, kw...)
    geom1 = first(Base.skipmissing(geoms))
    if GI.trait(geom1) isa GI.PointTrait
        if skipmissing
            Base.skipmissing((_extract_point_skip(A, g; kw...) for g in geoms))
        else
            (_extract_point(A, g; kw...) for g in geoms)
        end
    elseif GI.isgeometry(geom1) || GI.isfeature(geom1) || GI.isfeaturecollection(geom1)
        vals = (_extract(A, g; skipmissing, kw...) for g in geoms)
        skipmissing ? Base.skipmissing(vals) : vals
    else
        throw(ArgumentError("`data` does not contain geometry objects"))
    end
end
function _extract(A::RasterStackOrArray, ::GI.AbstractFeatureTrait, feature; kw...)
    _extract(A, GI.geometry(feature); kw...)
end
function _extract(A::RasterStackOrArray, ::GI.AbstractMultiPointTrait, geom; skipmissing=false, kw...)
    if skipmissing
        skipmissing(_extract_point_skip(A, p; kw...) for g in geoms)
    else
        (_extract_point(A, g; kw...) for g in geoms)
    end
end
function _extract(A::RasterStackOrArray, ::GI.AbstractGeometryTrait, geom; 
    names, geometry=true, index=false, skipmissing=false, kw...
)
    B = boolmask(geom; to=commondims(A, DEFAULT_POINT_ORDER), kw...)
    return _extract_geom(A, B, names; skipmissing, geometry, index)
end

@inline function _extract_geom(A, B, names; skipmissing, index, geometry)
    if skipmissing
        Base.skipmissing(_skipped_get(A, names, I; index, geometry) for I in CartesianIndices(B) if B[I])
    else
        (_normal_get(A, names, I; index, geometry) for I in CartesianIndices(B) if B[I])
    end
end

@inline _skipped_get(A, names, I; kw...) = _skipped_get(A, names, _prop_nt(A, I, names), I,; kw...)
@inline function _skipped_get(st::RasterStack, names, vals, I; geometry, index) 
    if any(x -> map((x, m) -> ismissing(x) || s === m, vals, missingval(st)), vals) 
        missing 
    else
        _fill_row(A, vals, I; geometry, index)
    end
end
@inline function _skipped_get(A, names, vals, I; geometry, index) 
    if any(x -> ismissing(x) || x === missingval(A), vals) 
        missing 
    else
        _fill_row(A, vals, I; geometry, index)
    end
end

@inline _normal_get(A, names, I; kw...) = _fill_row(A, _prop_nt(A, I, names), I; kw...)

_prop_nt(st::AbstractRasterStack, I, names::NamedTuple{K}) where K = t[I][K]
_prop_nt(A::AbstractRaster, I, names::NamedTuple{K}) where K = NamedTuple{K}((A[I],))

function _extract(x::RasterStackOrArray, ::GI.PointTrait, point; skipmissing=false)
    if skipmissing
        _extract_point(x, point; kw...)
    else
        _extract_point_skip(x, point; kw...)
    end
end
function _extract_point_skip(x::RasterStackOrArray, point; 
    dims, names, atol=nothing, geometry=true, index=false, kw...
)  
    ismissing(point) && return missing
    # Get the actual dimensions available in the object
    coords = map(DD.commondims(x, dims)) do d
        _dimcoord(d, point)
    end

    # Extract the values
    if any(map(ismissing, coords))
        return missing
    else
        selectors = map(dims, coords) do d, c
            _at_or_contains(d, c, atol)
        end
        if DD.hasselection(x, selectors)
            I = DD.dims2indices(x, selectors)
            val = x[I...]
            layer_vals = if x isa Raster
                NamedTuple{keys(names)}((val,))
            else
                val
            end
            return _fill_row(x, layer_vals, point, I; geometry, index)
        else
            return missing
        end
    end
end

function _extract_point(x::RasterStackOrArray, point; 
    dims, names::NamedTuple{K}, atol=nothing, geometry=true, index=false, kw...
) where K
    if ismissing(point) 
        layer_vals = map(_ -> missing, names)
        geom = missing
        I = missing
    else
        # Get the actual dimensions available in the object
        coords = map(DD.commondims(x, dims)) do d
            _dimcoord(d, point)
        end
        # Extract the values
        if any(map(ismissing, coords))
            # TODO test this branch somehow
            layer_vals = map(_ -> missing, names)
            geom = missing
            I = missing
        else
            selectors = map(dims, coords) do d, c
                _at_or_contains(d, c, atol)
            end
            if DD.hasselection(x, selectors)
                I = DD.dims2indices(x, selectors)
                layer_vals = x isa Raster ? NamedTuple{K}((x[I...],)) : x[I...][K]
            else
                I = missing
                layer_vals = map(_ -> missing, names)
            end
            geom = point
        end
    end

    return _fill_row(x, layer_vals, geom, I; geometry, index)
end

@inline function _fill_row(A, layer_vals::NamedTuple, I; kw...)
    _fill_row(A, layer_vals, DimPoints(A)[I], I; kw...)
end
@inline function _fill_row(A, layer_vals::NamedTuple, point, I; geometry, index)
    if geometry
        if index
            merge((; geometry=point, index=I), layer_vals)
        else
            merge((; geometry=point), layer_vals)
        end
    else
        if index
            merge((; index=I), layer_vals)
        else
            layer_vals
        end
    end
end

_names(A::AbstractRaster) = NamedTuple{(Symbol(name(A)),)}((Symbol(name(A)),))
_names(A::AbstractRasterStack) = NamedTuple{keys(A)}(keys(A))
