using DimensionalData.Lookups: _True, _False

_booltype(x) = x ? _True() : _False()
istrue(::_True) = true
istrue(::_False) = false

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

- `geometry`: include `:geometry` in retured `NamedTuple`, `true` by default.
- `index`: include `:index` of the `CartesianIndex` in retured `NamedTuple`, `false` by default.
- `name`: a `Symbol` or `Tuple` of `Symbol` corresponding to layer/s of a `RasterStack` to extract. All layers by default.
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
@inline function extract(x::RasterStackOrArray, data;
   name=_names(x), skipmissing=false, geometry=true, index=false, kw...
)
    n = DD._astuple(name)
    _extract(x::RasterStackOrArray, data, NamedTuple{n}(n); 
         # These keywords are converted to _True/_False for type stability later on
         # The @inline above helps constant propagation of the Bools
         geometry=_booltype(geometry), 
         index=_booltype(index), 
         skipmissing=_booltype(skipmissing), 
         kw...
    )
end
function _extract(x::RasterStackOrArray, data, names::NamedTuple{Names};
    dims=DD.dims(x, DEFAULT_POINT_ORDER), kw...
) where Names
    if !(data isa AbstractVector{<:GeoInterface.NamedTuplePoint}) && Tables.istable(data)
        geomcolnames = GI.geometrycolumns(data)
        if isnothing(geomcolnames) || !in(first(geomcolnames), Tables.columnnames(Tables.columns(data)))
            throw(ArgumentError("No `:geometry` column and `GeoInterface.geometrycolums(::$(typeof(data)))` does not define alternate columns"))
        end
        geometries = Tables.getcolumn(Tables.columns(data), first(geomcolnames))
        _extract(x, geometries; dims, names, kw...)
     else
        _extract(x, data; dims, names, kw...)
    end
end

function _extract(A::RasterStackOrArray, geom::Missing; names, kw...)
    T = _rowtype(A, geom; names, kw...)
    [_maybe_add_fields(T, A, map(_ -> missing, names), missing, missing)]
end
_extract(A::RasterStackOrArray, geom; kw...) = 
    _extract(A, GI.geomtrait(geom), geom; kw...)
_extract(A::RasterStackOrArray, ::Nothing, geom; kw...) =
    throw(ArgumentError("$geom is not a valid GeoInterface.jl geometry"))
function _extract(A::RasterStackOrArray, ::Nothing, geoms::AbstractArray; 
    names, skipmissing, kw...
)
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
        T = _rowtype(A, geom1; names, skipmissing, kw...)
        rows = T[]
        sizehint!(rows, length(geoms))
        for g in geoms
            e = _extract_point(T, A, g; names, skipmissing, kw...)
            ismissing(e) || push!(rows, e)
        end
    else
        # This will be a list of vectors, we need to flatten it into one
        rows = Iterators.flatten(_extract(A, g; names, skipmissing, kw...) for g in geoms)
        return istrue(skipmissing) ? collect(_skip_missing_rows(rows, _missingval_or_missing(A), names)) : collect(rows)
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
    T = _rowtype(A, GI.getpoint(geom, 1); names, skipmissing, kw...)
    rows = (_extract_point(T, A, p; skipmissing, kw...) for p in GI.getpoint(geom))
    return skipmissing isa _True ? collect(_skip_missing_rows(rows, _missingval_or_missing(A), names)) : collect(rows)
end
function _extract(A::RasterStackOrArray, t::GI.AbstractGeometryTrait, geom;
    names, skipmissing, kw...
)
    T = _rowtype(A, geom; names, skipmissing, kw...)
    # Make a raster mask of the geometry
    template = view(A, Touches(GI.extent(geom)))
    ods = otherdims(A, DEFAULT_POINT_ORDER)
    if length(ods) > 0
        template = view(template, map(d -> rebuild(d, firstindex(d)), ods)) 
    end
    B = boolmask(geom; to=template, kw...)
    offset = CartesianIndex(map(x -> first(x) - 1, parentindices(parent(template))))
    # Add a row for each pixel that is `true` in the mask
    rows = (_maybe_add_fields(T, A, _prop_nt(A, I + offset, names), I + offset) for I in CartesianIndices(B) if B[I])
    # Maybe skip missing rows
    return skipmissing ? collect(_skip_missing_rows(rows, _missingval_or_missing(A), names)) : collect(rows)
end
_extract(::Type{T}, x::RasterStackOrArray, trait::GI.PointTrait, point; kw...) where T =
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


@inline _prop_nt(st::AbstractRasterStack, I, ::NamedTuple{K}) where K = st[I][K]
@inline _prop_nt(A::AbstractRaster, I, ::NamedTuple{K}) where K = NamedTuple{K}((A[I],))

# Extract a single point
Base.@assume_effects :foldable function _extract_point(::Type{T}, x::RasterStackOrArray, point::Missing;
    names::NamedTuple{K}, skipmissing, kw...
) where {T,K}
    if istrue(skipmissing)
        missing
    else
        layer_vals = map(_ -> missing, names)
        geom = missing
        I = missing
        return _maybe_add_fields(T, x, layer_vals, geom, I)
    end
end
Base.@assume_effects :foldable function _extract_point(::Type{T}, x::RasterStackOrArray, point::Tuple;
    dims, names::NamedTuple{K}, atol=nothing, kw...
) where {T,K}
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
        inds = DD.selectindices(x, selectors; err=_False())
        if isnothing(inds)
            # Otherwise return `missing` for everything
            I = missing
            layer_vals = map(_ -> missing, names)
        else
            I = inds
            layer_vals = x isa Raster ? NamedTuple{K}((x[I...],)) : x[I...][K]
        end
        # Return the point for the geometry, either way
        geom = point
    end

    return _maybe_add_fields(T, x, layer_vals, geom, I)
end

# Maybe add optional fields
@inline function _maybe_add_fields(::Type{T}, A, layer_vals::NamedTuple, I) where T
    _maybe_add_fields(T, A, layer_vals, DimPoints(A)[I], I)::T
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

@inline _nametypes(::Raster{T}, ::NamedTuple{Names}, skipmissing::_True) where {T,Names} = (T,)
@inline _nametypes(::Raster{T}, ::NamedTuple{Names}, skipmissing::_False) where {T,Names} = (Union{Missing,T},)
# @inline _nametypes(::RasterStack{K,T}, ::NamedTuple{Names}) = NamedTuple{T}()

# _rowtype returns the complete NamedTuple type for a point row
# This code is entrirely for types stability and performance.
_rowtype(x, g; geometry, index, skipmissing, names, kw...) = _rowtype(x, g, geometry, index, skipmissing, names)
function _rowtype(x, g::G, geometry::_True, index::_True, skipmissing::_True, names::NamedTuple{Names}) where {G,Names}
    keys = (:geometry, :index, Names...,) 
    types = Tuple{typeof(g),Tuple{Int,Int},_nametypes(x, names, skipmissing)...}
    NamedTuple{keys,types}
end
function _rowtype(x, g::G, geometry::_True, index::_False, skipmissing::_True, names::NamedTuple{Names}) where {G,Names}
    keys = (:geometry, Names...,)
    types = Tuple{typeof(g),_nametypes(x, names, skipmissing)...}
    NamedTuple{keys,types}
end
function _rowtype(x, g::G, geometry::_False, index::_True, skipmissing::_True, names::NamedTuple{Names}) where {G,Names}
    keys = (:index, Names...,) 
    types = Tuple{Tuple{Int,Int},_nametypes(x, names, skipmissing)...}
    NamedTuple{keys,types}
end
function _rowtype(x, g::G, geometry::_False, index::_False, skipmissing::_True, names::NamedTuple{Names}) where {G,Names}
    keys = Names
    types = Tuple{_nametypes(x, names, skipmissing)...}
    NamedTuple{keys,types}
end
function _rowtype(x, g::G, geometry::_True, index::_True, skipmissing::_False, names::NamedTuple{Names}) where {G,Names}
    keys = (:geometry, :index, names...,) 
    types = Tuple{Union{Missing,typeof(g)},Union{Missing,Tuple{Int,Int}},_nametypes(x, names)...}
    NamedTuple{keys,types}
end
function _rowtype(x, g::G, geometry::_True, index::_False, skipmissing::_False, names::NamedTuple{Names}) where {G,Names}
    keys = (:geometry, Names...,)
    types = Tuple{Union{Missing,typeof(g)},_nametypes(x, names, skipmissing)...}
    NamedTuple{keys,types}
end
function _rowtype(x, g::G, geometry::_False, index::_True, skipmissing::_False, names::NamedTuple{Names}) where {G,Names}
    keys = (:index, Names...,) 
    types = Tuple{Union{Missing,Tuple{Int,Int}},_nametypes(x, names, skipmissing)...}
    NamedTuple{keys,types}
end
function _rowtype(x, g::G, geometry::_False, index::_False, skipmissing::_False, names::NamedTuple{Names}) where {G,Names}
    keys = Names
    types = Tuple{_nametypes(x, names, skipmissing)...}
    NamedTuple{keys,types}
end
