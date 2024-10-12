"""
    extract(x, data; kw...)

Extracts the value of `Raster` or `RasterStack` at given points, returning
an iterable of `NamedTuple` with properties for `:geometry` and raster or
stack layer values.

Note that if objects have more dimensions than the length of the point tuples,
sliced arrays or stacks will be returned instead of single values.

# Arguments

- `x`: a `Raster` or `RasterStack` to extract values from.
$DATA_ARGUMENT

# Keywords

- `geometry`: include `:geometry` in returned `NamedTuple`, `true` by default.
- `index`: include `:index` of the `CartesianIndex` in returned `NamedTuple`, `false` by default.
- `name`: a `Symbol` or `Tuple` of `Symbol` corresponding to layer/s of a `RasterStack` to extract. All layers by default.
- `skipmissing`: skip missing points automatically.
- `atol`: a tolerance for floating point lookup values for when the `Lookup`
    contains `Points`. `atol` is ignored for `Intervals`.
-$GEOMETRYCOLUMN_KEYWORD

# Example

Here we extract points matching the occurrence of the Mountain Pygmy Possum,
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
@inline function extract(x::RasterStackOrArray, data;
   names=_names(x), name=names, skipmissing=false, geometry=true, index=false, geometrycolumn=nothing, kw...
)
    n = DD._astuple(name)
    _extract(x, data;
         dims=DD.dims(x, DEFAULT_POINT_ORDER),
         names=NamedTuple{n}(n),
         # These keywords are converted to _True/_False for type stability later on
         # The @inline above helps constant propagation of the Bools
         geometry=_booltype(geometry), 
         index=_booltype(index), 
         skipmissing=_booltype(skipmissing), 
         geometrycolumn,
         kw...
    )
end

function _extract(A::RasterStackOrArray, geom::Missing, names, kw...)
    T = _extractrowtype(A, geom; names, kw...)
    [_maybe_add_fields(T, map(_ -> missing, names), missing, missing)]
end
function _extract(A::RasterStackOrArray, geom; names, kw...)
    _extract(A, GI.geomtrait(geom), geom; names, kw...)
end
function _extract(A::RasterStackOrArray, ::Nothing, data; 
    names, skipmissing, geometrycolumn, kw...
)
    geoms = _get_geometries(data, geometrycolumn)
    T = if istrue(skipmissing)
        _extractrowtype(A, nonmissingtype(eltype(geoms)); names, skipmissing, kw...)
    else
        _extractrowtype(A, eltype(geoms); names, skipmissing, kw...)
    end
    # Handle empty / all missing cases
    (length(geoms) > 0 && any(!ismissing, geoms)) || return T[]
    
    geom1 = first(Base.skipmissing(geoms))
    trait1 = GI.trait(geom1)
    # We need to split out points from other geoms
    # TODO this will fail with mixed point/geom vectors
    if trait1 isa GI.PointTrait
        rows = Vector{T}(undef, length(geoms))
        if istrue(skipmissing)
            j = 1
            for i in eachindex(geoms)
                g = geoms[i]
                ismissing(g) && continue
                e = _extract_point(T, A, g, skipmissing; names, kw...)
                if !ismissing(e) 
                    rows[j] = e
                    j += 1
                end
                nothing
            end
            deleteat!(rows, j:length(rows))
        else
            for i in eachindex(geoms)
                g = geoms[i]
                rows[i] = _extract_point(T, A, g, skipmissing; names, kw...)::T
                nothing
            end
        end
        return rows
    else
        # This will be a list of vectors, we need to flatten it into one
        rows = Iterators.flatten(_extract(A, g; names, skipmissing, kw...) for g in geoms)
        if istrue(skipmissing) 
            return collect(_skip_missing_rows(rows, _missingval_or_missing(A), names)) 
        else
            return collect(rows)
        end
    end
end
function _extract(A::RasterStackOrArray, ::GI.AbstractFeatureTrait, feature; kw...)
    _extract(A, GI.geometry(feature); kw...)
end
function _extract(A::RasterStackOrArray, ::GI.FeatureCollectionTrait, fc; kw...)
    # Fall back to the Array/iterator method for feature collections
    _extract(A, [GI.geometry(f) for f in GI.getfeature(fc)]; kw...)
end
function _extract(A::RasterStackOrArray, ::GI.AbstractMultiPointTrait, geom; 
    skipmissing, kw...
)
    T = _extractrowtype(A, GI.getpoint(geom, 1); names, skipmissing, kw...)
    rows = (_extract_point(T, A, p, skipmissing; kw...) for p in GI.getpoint(geom))
    return skipmissing isa _True ? collect(_skip_missing_rows(rows, _missingval_or_missing(A), names)) : collect(rows)
end
function _extract(A::RasterStackOrArray, ::GI.PointTrait, geom; 
    skipmissing, kw...
)
    T = _extractrowtype(A, geom; names, skipmissing, kw...)
    _extract_point(T, A, geom, skipmissing; kw...)
end
function _extract(A::RasterStackOrArray, t::GI.AbstractGeometryTrait, geom;
    names, skipmissing, kw...
)
    # Make a raster mask of the geometry
    template = view(A, Touches(GI.extent(geom)))
    ods = otherdims(A, DEFAULT_POINT_ORDER)
    if length(ods) > 0
        template = view(template, map(d -> rebuild(d, firstindex(d)), ods)) 
    end
    p = first(GI.getpoint(geom))
    tuplepoint = if GI.is3d(geom) 
        GI.x(p), GI.y(p), GI.z(p)
    else
        GI.x(p), GI.y(p)
    end
    T = _extractrowtype(A, tuplepoint; names, skipmissing, kw...)
    B = boolmask(geom; to=template, kw...)
    offset = CartesianIndex(map(x -> first(x) - 1, parentindices(parent(template))))
    # Add a row for each pixel that is `true` in the mask
    rows = (_missing_or_fields(T, A, Tuple(I + offset), names, skipmissing) for I in CartesianIndices(B) if B[I])
    # Maybe skip missing rows
    if istrue(skipmissing)
        return collect(Base.skipmissing(rows)) 
    else
        return collect(rows)
    end
end
_extract(::Type{T}, A::RasterStackOrArray, trait::GI.PointTrait, point; skipmissing, kw...) where T =
    _extract_point(T, A, point, skipmissing; kw...)

function _missing_or_fields(::Type{T}, A, I, names, skipmissing) where T
    props = _prop_nt(A, I, names)
    if istrue(skipmissing) && _ismissingval(A, props) 
        missing
    else
        geom = DimPoints(A)[I...]
        _maybe_add_fields(T, props, geom, I)
    end
end

_ismissingval(A::Union{Raster,RasterStack}, props) = 
    _ismissingval(missingval(A), props)
_ismissingval(A::Union{Raster,RasterStack}, props::NamedTuple) = 
    _ismissingval(missingval(A), props)
_ismissingval(mvs::NamedTuple, props::NamedTuple{K}) where K = 
    any(k -> ismissing(props[k]) || props[k] === mvs[k], K)
_ismissingval(mv, props::NamedTuple) = any(x -> ismissing(x) || x === mv, props)
_ismissingval(mv, prop) = (mv === prop)

@inline _prop_nt(st::AbstractRasterStack, I, ::NamedTuple{K}) where K = st[I...][K]
@inline _prop_nt(A::AbstractRaster, I, ::NamedTuple{K}) where K = NamedTuple{K}((A[I...],))

# Extract a single point
@inline function _extract_point(::Type, x::RasterStackOrArray, point::Missing, skipmissing::_True; kw...)
    return missing
end
@inline function _extract_point(::Type{T}, x::RasterStackOrArray, point::Missing, skipmissing::_False;
    names, kw...
) where T
    props = map(_ -> missing, names)
    geom = missing
    I = missing
    return _maybe_add_fields(T, props, geom, I)
end
@inline function _extract_point(::Type{T}, x::RasterStackOrArray, point, skipmissing::_True;
    dims, names::NamedTuple{K}, atol=nothing, kw...
) where {T,K}
    # Get the actual dimensions available in the object
    coords = map(DD.commondims(x, dims)) do d
        _dimcoord(d, point)
    end
    # If any are coordinates missing, also return missing for everything
    if any(map(ismissing, coords))
        return missing
    else
        selector_dims = map(dims, coords) do d, c
            _at_or_contains(d, c, atol)
        end
        selectors = map(val, DD.dims(selector_dims, DD.dims(x)))
        # Check the selector is in bounds / actually selectable
        I = DD.selectindices(DD.dims(x, dims), selectors; err=_False())::Union{Nothing,Tuple{Int,Int}}
        if isnothing(I)
            return missing
        else
            D = map(rebuild, selector_dims, I)
            if x isa Raster 
                prop = x[D]
                _ismissingval(missingval(x), prop) && return missing
                props = NamedTuple{K,Tuple{eltype(x)}}((prop,))
            else
                props = x[D][K]
                _ismissingval(missingval(x), props) && return missing
            end
            return _maybe_add_fields(T, props, point, I)
        end
    end
end
@inline function _extract_point(::Type{T}, x::RasterStackOrArray, point, skipmissing::_False;
    dims, names::NamedTuple{K}, atol=nothing, kw...
) where {T,K}
    # Get the actual dimensions available in the object
    coords = map(DD.commondims(x, dims)) do d
        _dimcoord(d, point)
    end
    # If any are coordinates missing, also return missing for everything
    if any(map(ismissing, coords))
        return _maybe_add_fields(T, map(_ -> missing, names), missing, missing)
    else
        selector_dims = map(dims, coords) do d, c
            _at_or_contains(d, c, atol)
        end
        selectors = map(val, DD.dims(selector_dims, DD.dims(x)))
        # Check the selector is in bounds / actually selectable
        I = DD.selectindices(DD.dims(x, dims), selectors; err=_False())::Union{Nothing,Tuple{Int,Int}}
        if isnothing(I)
            return _maybe_add_fields(T, map(_ -> missing, names), point, missing)
        else
            D = map(rebuild, selector_dims, I)
            props = if x isa Raster 
                NamedTuple{K,Tuple{eltype(x)}}((x[D],))
            else
                NamedTuple(x[D])[K]
            end
            return _maybe_add_fields(T, props, point, I)
        end
    end
end

# Maybe add optional fields
Base.@assume_effects :total function _maybe_add_fields(::Type{T}, props::NamedTuple, point, I)::T where {T<:NamedTuple{K}} where K
    if :geometry in K
        :index in K ? merge((; geometry=point, index=I), props) : merge((; geometry=point), props)
    else
        :index in K ? merge((; index=I), props) : props
    end
end

function _extractrowtype(x, g; geometry, index, skipmissing, names, kw...)
    I = if istrue(skipmissing)
        Tuple{Int, Int}
    else
        Union{Missing, Tuple{Int, Int}}
    end
    _extractrowtype(x, g, I; geometry, index, skipmissing, names, kw...)
end
function _extractrowtype(x, g, ::Type{I}; geometry, index, skipmissing, names, kw...) where I
    G = if istrue(skipmissing)
        nonmissingtype(typeof(g))
    else
        typeof(g)
    end
    _extractrowtype(x, G, I; geometry, index, skipmissing, names, kw...)
end
_extractrowtype(x, ::Type{G}, ::Type{I}; geometry, index, skipmissing, names, kw...) where {G, I} =
    _rowtype(x, G, I; geometry, index, skipmissing, names)

@inline _skip_missing_rows(rows, ::Missing, names) = 
    Iterators.filter(row -> !any(ismissing, row), rows)
@inline _skip_missing_rows(rows, missingval, names) = 
    Iterators.filter(row -> !any(x -> ismissing(x) || x === missingval, row), rows)
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
