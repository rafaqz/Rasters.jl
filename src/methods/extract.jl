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
- `index`: include an `:index` column of the `CartesianIndex` for each value, `false` by default.
- `names`: `Tuple` of `Symbol` corresponding to layers of a `RasterStack`. All layers by default.
- `skipmissing`: skip missing points automatically.
- `atol`: a tolorerance for floating point lookup values for when the `Lookup`
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
extract(st, obs; skipmissing=true)

# output
5-element Vector{NamedTuple{(:geometry, :bio1, :bio3, :bio5, :bio7, :bio12)}}:
 (geometry = (0.21, 40.07), bio1 = 17.077084f0, bio3 = 41.20417f0, bio5 = 30.1f0, bio7 = 24.775f0, bio12 = 446.0f0)
 (geometry = (0.03, 39.97), bio1 = 17.076923f0, bio3 = 39.7983f0, bio5 = 29.638462f0, bio7 = 24.153847f0, bio12 = 441.0f0)
 (geometry = (0.03, 39.97), bio1 = 17.076923f0, bio3 = 39.7983f0, bio5 = 29.638462f0, bio7 = 24.153847f0, bio12 = 441.0f0)
 (geometry = (0.52, 40.37), bio1 = missing, bio3 = missing, bio5 = missing, bio7 = missing, bio12 = missing)
 (geometry = (0.32, 40.24), bio1 = 16.321388f0, bio3 = 41.659454f0, bio5 = 30.029825f0, bio7 = 25.544561f0, bio12 = 480.0f0)
```

Note: passing in arrays, geometry collections or feature collections 
containing a mix of points and other geometries has undefined results.
"""
function extract end
function extract(x::RasterStackOrArray, data;
    dims=DD.dims(x, DEFAULT_POINT_ORDER), names=_names(x), geometry=true, index=false, kw...
)
    T = if geometry 
        keys = index ? (:geometry, :index, names...,) : (:geometry, names...,)
        NamedTuple{keys}
    else
        keys = index ? (:index, names...,) : (names...,)
        NamedTuple{keys}
    end
    names = NamedTuple{names}(names)
    if !(data isa AbstractVector{<:GeoInterface.NamedTuplePoint}) && Tables.istable(data)
        geomcolnames = GI.geometrycolumns(data)
        if isnothing(geomcolnames) || !in(first(geomcolnames), Tables.columnnames(Tables.columns(data)))
            throw(ArgumentError("No `:geometry` column and `GeoInterface.geometrycolums(::$(typeof(data)))` does not define alternate columns"))
        end
        geometries = Tables.getcolumn(Tables.columns(data), first(geomcolnames))
        _extract(T, x, geometries; dims, names, kw...)
     else
        _extract(T, x, data; dims, names, kw...)
    end
end

function _extract(T, A::RasterStackOrArray, geom::Missing; names, kw...)
    [_maybe_add_fields(T, A, map(_ -> missing, names), missing, missing)]
end
_extract(T, A::RasterStackOrArray, geom; kw...) = 
    _extract(T, A, GI.geomtrait(geom), geom; kw...)
_extract(T, A::RasterStackOrArray, ::Nothing, geom; kw...) =
    throw(ArgumentError("$geom is not a valid GeoInterface.jl geometry"))
function _extract(T, A::RasterStackOrArray, ::Nothing, geoms::AbstractArray; names, skipmissing=false, kw...)
    # Handle empty / all missing cases
    (length(geoms) > 0 && any(!ismissing, geoms)) || return T[]

    # Handle cases with some invalid geometries
    invalid_geom_idx = findfirst(g -> !ismissing(g) && GI.geomtrait(g) === nothing, geoms)
    invalid_geom_idx === nothing || throw(ArgumentError("$(geoms[invalid_geom_idx]) is not a valid GeoInterface.jl geometry"))

    geom1 = first(Base.skipmissing(geoms))
    trait1 = GI.trait(geom1)
    # We need to split out points from other geoms
    # TODO this will fail with mixed point/geom vectors
    if trait1 isa GI.PointTrait
        rows = (_extract_point(T, A, g; names, kw...) for g in geoms)
        return skipmissing ? collect(_skip_missing_rows(rows, _missingval_or_missing(A), names)) : collect(rows)
    else
        # This will be a list of vectors, we need to flatten it into one
        rows = Iterators.flatten(_extract(T, A, g; names, skipmissing, kw...) for g in geoms)
        return skipmissing ? collect(_skip_missing_rows(rows, _missingval_or_missing(A), names)) : collect(rows)
    end
end
function _extract(T, A::RasterStackOrArray, ::GI.AbstractFeatureTrait, feature; kw...)
    _extract(T, A, GI.geometry(feature); kw...)
end
function _extract(T, A::RasterStackOrArray, ::GI.FeatureCollectionTrait, fc; kw...)
    # Fall back to the Array/iterator method for feature collections
    _extract(T, A, [GI.geometry(f) for f in GI.getfeature(fc)]; kw...)
end
function _extract(T, A::RasterStackOrArray, ::GI.AbstractMultiPointTrait, geom; skipmissing=false, kw...)
    rows = (_extract_point(T, A, g; kw...) for g in GI.getpoint(geom))
    return skipmissing ? collect(_skip_missing_rows(rows, _missingval_or_missing(A), names)) : collect(rows)
end
function _extract(T, A::RasterStackOrArray, ::GI.AbstractGeometryTrait, geom;
    names, skipmissing=false, kw...
)
    # Make a raster mask of the geometry
    B = boolmask(geom; to=commondims(A, DEFAULT_POINT_ORDER), kw...)
    # Add a row for each pixel that is `true` in the mask
    rows = (_maybe_add_fields(T, A, _prop_nt(A, I, names), I) for I in CartesianIndices(B) if B[I])
    # Maybe skip missing rows
    return skipmissing ? collect(_skip_missing_rows(rows, _missingval_or_missing(A), names)) : collect(rows)
end
_extract(T, x::RasterStackOrArray, trait::GI.PointTrait, point; kw...) =
    _extract_point(T, x, point; kw...)

@inline _skip_missing_rows(rows, ::Missing, names) = Iterators.filter(row -> !any(ismissing, row), rows)
@inline _skip_missing_rows(rows, missingval, names) = Iterators.filter(row -> !any(x -> ismissing(x) || x === missingval, row), rows)
@inline function _skip_missing_rows(rows, missingval::NamedTuple, names::NamedTuple{K}) where K
    # first check if all fields are equal - if so just call with the first value
    if Base.allequal(missingval) == 1
        return _skip_missing_rows(rows, first(missingval), names)
    else
        Iterators.filter(rows) do row
            # rows may or may not contain a :geometry field, so map over keys instead
            !any(key -> ismissing(row[key]) || row[key] === missingval[key], K)
        end
    end
end


@inline _prop_nt(st::AbstractRasterStack, I, names::NamedTuple{K}) where K = st[I][K]
@inline _prop_nt(A::AbstractRaster, I, names::NamedTuple{K}) where K = NamedTuple{K}((A[I],))

# Extract a single point
function _extract_point(T, x::RasterStackOrArray, point;
    dims, names::NamedTuple{K}, atol=nothing, kw...
) where K
    # The point itself might be missing, so return missing for every field
    if ismissing(point)
        layer_vals = map(_ -> missing, names)
        geom = missing
        I = missing
    else
        # Get the actual dimensions available in the object
        coords = map(DD.commondims(x, dims)) do d
            _dimcoord(d, point)
        end
        # If any are coordinates missing, also return missing for everything
        if any(map(ismissing, coords))
            layer_vals = map(_ -> missing, names)
            geom = missing
            I = missing
        else
            selectors = map(dims, coords) do d, c
                _at_or_contains(d, c, atol)
            end
            # Check the selector is in bounds / actually selectable
            if DD.hasselection(x, selectors)
                I = DD.dims2indices(x, selectors)
                layer_vals = x isa Raster ? NamedTuple{K}((x[I...],)) : x[I...][K]
            else
                # Otherwise return `missing` for everything
                I = missing
                layer_vals = map(_ -> missing, names)
            end
            # Return the point for the geometry, either way
            geom = point
        end
    end

    return _maybe_add_fields(T, x, layer_vals, geom, I)
end
function _extract_point(T, A::RasterStackOrArray, point::Missing; names, kw...)
    # Missing points return a single row
    return _maybe_add_fields(T, A, map(_ -> missing, names), missing, missing)
end

# Maybe add optional fields
@inline function _maybe_add_fields(T, A, layer_vals::NamedTuple, I)
    _maybe_add_fields(T, A, layer_vals, DimPoints(A)[I], I)
end
@inline function _maybe_add_fields(::Type{T}, A, layer_vals::NamedTuple, point, I)::T where {T<:NamedTuple{K}} where K
    if :geometry in K
        :index in K ? merge((; geometry=point, index=I), layer_vals) : merge((; geometry=point), layer_vals)
    else
        :index in K ? merge((; index=I), layer_vals) : layer_vals
    end
end

_names(A::AbstractRaster) = (Symbol(name(A)),)
_names(A::AbstractRasterStack) = keys(A)
