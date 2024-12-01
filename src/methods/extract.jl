# Object to hold extract keywords
struct Extractor{T,P,K,N<:NamedTuple{K},G,I,S,F}
    names::N
    geometry::G
    index::I
    skipmissing::S
    flatten::F
end
function Extractor(x, data;
    name::NTuple{<:Any,Symbol},
    skipmissing,
    flatten,
    geometry,
    index,
    kw...
)
    nt = NamedTuple{name}(name)
    g = _booltype(geometry)
    i = _booltype(index)
    sm = _booltype(skipmissing)
    f = _booltype(flatten)
    T = _geom_rowtype(x, data; geometry=g, index=i, names=nt, skipmissing=sm)
    P = _proptype(x; skipmissing=sm, names=nt)
    Extractor{T,P,name,typeof(nt),typeof(g),typeof(i),typeof(sm),typeof(f)}(nt, g, i, sm, f)
end

Base.eltype(::Extractor{T}) where T = T

mutable struct LineRefs{T}
    i::Int
    prev::CartesianIndex{2}
    rows::Vector{T}
    function LineRefs{T}() where T
        i = 1
        prev = CartesianIndex((typemin(Int), typemin(Int)))
        new{T}(i, prev, Vector{T}())
    end
end

function _initialise!(lr::LineRefs)
    lr.i = 1
    lr.prev = CartesianIndex((typemin(Int), typemin(Int)))
end

_geom_rowtype(A, geom; kw...) =
    _geom_rowtype(A, GI.trait(geom), geom; kw...)
_geom_rowtype(A, ::Nothing, geoms; kw...) =
    _geom_rowtype(A, first(geoms); kw...)
_geom_rowtype(A, ::GI.AbstractGeometryTrait, geom; kw...) =
    _geom_rowtype(A, first(GI.getpoint(geom)); kw...)
function _geom_rowtype(A, ::GI.AbstractPointTrait, p; kw...)
    tuplepoint = if GI.is3d(p)
        (GI.x(p), GI.y(p), GI.z(p))
    else
        (GI.x(p), GI.y(p))
    end
    return _rowtype(A, tuplepoint; kw...)
end

# skipinvalid: can G and I be missing. skipmissing: can nametypes be missing
function _proptype(x;
    skipmissing, names::NamedTuple{K}, kw...
) where K
    NamedTuple{K,Tuple{_nametypes(x, names, skipmissing)...}}
end

"""
    extract(x, geometries; kw...)

Extracts the value of `Raster` or `RasterStack` for the passed in geometries, 
returning an `Vector{NamedTuple}` with properties for `:geometry` and `Raster`
or `RasterStack` layer values.

For lines, linestrings and linear rings points are extracted for each pixel that
the line touches. 

For polygons, all cells witih centers covered by the polygon are returned.

Note that if objects have more dimensions than the length of the point tuples,
sliced arrays or stacks will be returned instead of single values.

# Arguments

- `x`: a `Raster` or `RasterStack` to extract values from.
$DATA_ARGUMENT

# Keywords

- `geometry`: include `:geometry` column in rows. `true` by default.
- `index`: include `:index` of the `CartesianIndex` in rows, `false` by default.
- `name`: a `Symbol` or `Tuple` of `Symbol` corresponding to layer/s of a `RasterStack` to extract.
    All layers are extracted by default.
- `skipmissing`: skip missing points automatically.
- `flatten`: flatten extracted points from multiple geometries into a single 
    vector. `true` by default. Unmixed point geometries are always flattened.
- `atol`: a tolerance for floating point lookup values for when the `Lookup`
    contains `Points`. `atol` is ignored for `Intervals`.
$BOUNDARY_KEYWORD
$GEOMETRYCOLUMN_KEYWORD

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
    geometrycolumn=nothing,
    names::NTuple{<:Any,Symbol}=_names(x),
    name::NTuple{<:Any,Symbol}=names,
    skipmissing=false,
    flatten=true,
    geometry=true,
    index=false,
    kw...
)
    n = DD._astuple(name)
    g, g1 = if GI.isgeometry(data)
        data, data
    elseif GI.isfeature(data)
        g = GI.geometry(data)
        g, g
    else
        gs = _get_geometries(data, geometrycolumn)
        gs, first(Base.skipmissing(gs))
    end

    xp = _prepare_for_burning(x)
    e = Extractor(xp, g1; name, skipmissing, flatten, geometry, index)
    return _extract(xp, g, e; kw...)
end

# TODO use a GeometryOpsCore method like `applyreduce` here?
function _extract(A::RasterStackOrArray, geom::Missing, e; kw...)
    return if istrue(e.skipmissing)
        T[]
    else
        T[_maybe_add_fields(e, map(_ -> missing, e.names), missing, missing)]
    end
end
function _extract(A::RasterStackOrArray, geom, e; kw...)
    _extract(A, GI.geomtrait(geom), geom, e; kw...)
end
function _extract(A::RasterStackOrArray, ::Nothing, geoms::AbstractArray, e::Extractor{T};
    threaded=false, progress=true, kw...
) where T
    # Handle empty / all missing cases
    (length(geoms) > 0 && any(!ismissing, geoms)) || return T[]

    geom1 = first(skipmissing(geoms))
    trait1 = GI.trait(geom1)
    # We split out points from other geoms for performance
    if trait1 isa GI.PointTrait
        allpoints = true
        i = 1
        rows = Vector{T}(undef, length(geoms))
        if istrue(e.skipmissing)
            i = 1
            for i in eachindex(geoms)
                g = geoms[i]
                ismissing(g) && continue
                geomtrait(g) isa GI.PointTrait || break
                i += _extract_point!(rows, A, g, e, i; kw...)
            end
            deleteat!(rows, i:length(rows))
        else
            _run(1:length(geoms)) do i
                g = geoms[i]
                geomtrait(g) isa GI.PointTrait || return nothing
                _extract_point!(rows, A, g, e, i; kw...)::T
                return nothing
            end
        end
        # If we found a non-point geometry,
        # ignore these rows and start again generically
        allpoints && return rows
    end

    row_vecs = Vector{Vector{T}}(undef, length(geoms))
    _run(1:length(geoms), threaded, progress, "Extracting geometries...") do i
        row_vecs[i] = _extract(A, geoms[i], e; kw...)
    end
    if istrue(e.flatten)
        n = sum(map(length, row_vecs))
        rows = Vector{T}(undef, n)
        i = 1
        for rs in row_vecs
            for r in rs
                rows[i] = r 
                i += 1
            end
        end
        return rows
    else
        return row_vecs
    end
end
function _extract(A::RasterStackOrArray, ::GI.AbstractMultiPointTrait, geom, e; kw...)
    n = GI.npoint(geom)
    rows = _init_rows(e, n)
    for p in GI.getpoint(geom)
        i += _extract_point!(rows, A, e, p, i; kw...)
    end
    if istrue(skipmissing)
        # Remove excees rows where missing
        deleteat!(rows, i:length(rows))
    end
    return rows
end
function _extract(A::RasterStackOrArray, ::GI.PointTrait, p, e; kw...)
    _extract_point!(rows, A, e, p, length(rows); kw...)
end
@noinline function _extract(
    A::RasterStackOrArray, trait::GI.AbstractLineStringTrait, geom, e::Extractor{T};
    line_refs=LineRefs{T}(), 
    kw...
) where T
    _initialise!(line_refs)
    offset = CartesianIndex((0, 0))
    dp = DimPoints(A)
    function _length_callback(n) 
        rows = line_refs.rows
        resize!(rows, n + line_refs.i - 1)
    end

    _burn_lines!(_length_callback, dims(A), geom) do D
        I = CartesianIndex(map(val, D))
        # Avoid duplicates from adjacent line segments
        line_refs.prev == I && return nothing
        line_refs.prev = I
        # Make sure we only hit this pixel once
        # D is always inbounds
        line_refs.i += _maybe_set_row!(line_refs.rows, e.skipmissing, e, A, dp, offset, I, line_refs.i)
        return nothing
    end
    rows = line_refs.rows
    deleteat!(rows, line_refs.i+1:length(rows))
    return rows
end
@noinline function _extract(
    A::RasterStackOrArray, t::GI.AbstractGeometryTrait, geom, e; kw...
)
    # Make a raster mask of the geometry
    dims, offset = _template(A, geom)
    B = boolmask(geom; to=dims, burncheck_metadata=NoMetadata(), kw...)
    n = count(B)
    dp = DimPoints(A)
    i = 1
    rows = _init_rows(e, n)
    # Add a row for each pixel that is `true` in the mask
    for I in CartesianIndices(B)
        B[I] || continue
        i += _maybe_set_row!(rows, e.skipmissing, e, A, dp, offset, I, i)
    end
    # Cleanup
    deleteat!(rows, i:length(rows))
    return rows
end

function _template(x, geom)
    ods = otherdims(x, DEFAULT_POINT_ORDER)
    # Build a dummy raster in case this is a stack
    # views are just easier on an array than dims
    t1 = Raster(FillArrays.Zeros(size(x)), dims(x))
    t2 = if length(ods) > 0
        view(t1, map(d -> rebuild(d, firstindex(d)), ods))
    else
        t1
    end
    I = dims2indices(dims(t2), Touches(GI.extent(geom)))
    t3 = view(t2, I...)
    offset = CartesianIndex(map(first, I))
    return dims(t3), offset
end

Base.@assume_effects :foldable function _maybe_set_row!(
    rows, skipmissing::_True, e, A, dp, offset, I, i;
    props=_prop_nt(e, A, I)
)
    return if !_ismissingval(A, props)
        _maybe_set_row!(rows, _False(), e, A, dp, offset, I; props)
    else
        0
    end
end
Base.@assume_effects :foldable function _maybe_set_row!(
    rows::Vector{T}, skipmissing::_False, e, A, dp, offset, I, i;
    props=_prop_nt(e, A, I)
) where T
    Ioff = I + offset
    geom = dp[Ioff]
    i <= length(rows) || @show i length(rows)
    rows[i] = _maybe_add_fields(T, props, geom, Tuple(Ioff))
    return 1
end

_init_rows(e::Extractor{T}, n=0) where T = Vector{T}(undef, n)

Base.@assume_effects :foldable _ismissingval(A::Union{Raster,RasterStack}, props)::Bool =
    _ismissingval(missingval(A), props)
Base.@assume_effects :foldable _ismissingval(A::Union{Raster,RasterStack}, props::NamedTuple)::Bool =
    _ismissingval(missingval(A), props)
Base.@assume_effects :foldable _ismissingval(mvs::NamedTuple, props::NamedTuple)::Bool =
    any(DD.unrolled_map((p, mv) -> ismissing(p) || p === mv, props, mvs))
Base.@assume_effects :foldable _ismissingval(mv, props::NamedTuple)::Bool =
    any(DD.unrolled_map(p -> ismissing(p) || p === mv, props))
Base.@assume_effects :foldable _ismissingval(mv, prop)::Bool = (mv === prop)

# We always return NamedTuple
Base.@assume_effects :foldable _prop_nt(::Extractor{<:Any,P,K}, st::AbstractRasterStack, I) where {P,K} =
    P(st[K][I])::P
Base.@assume_effects :foldable _prop_nt(::Extractor{<:Any,P,K}, st::AbstractRasterStack{K}, I) where {P,K} =
    P(st[I])::P
Base.@assume_effects :foldable _prop_nt(::Extractor{<:Any,P,K}, A::AbstractRaster, I) where {P,K} =
    P((A[I],))::P

# Extract a single point
# Missing point to remove
@inline function _extract_point!(
    rows::Vector{T}, x::RasterStackOrArray, skipmissing::_True, point::Missing, i;
    kw...
) where T
    return 0
end
# Missing point to keep
@inline function _extract_point!(
    rows::Vector{T}, x::RasterStackOrArray, skipmissing::_False, point::Missing, i;
    names, kw...
) where T
    props = map(_ -> missing, names)
    geom = missing
    I = missing
    return _maybe_add_fields(T, props, geom, I)
end
# Normal point with missing / out of bounds data removed
@inline function _extract_point!(
    rows::Vector{T}, x::RasterStackOrArray, skipmissing::_True, point, i;
    dims=DD.dims(x, DEFAULT_POINT_ORDER),
    names::NamedTuple{K},
    atol=nothing,
    kw...
) where {T,K}
    # Get the actual dimensions available in the object
    coords = map(DD.commondims(x, dims)) do d
        _dimcoord(d, point)
    end
    # If any are coordinates missing, also return missing for everything
    if !any(map(ismissing, coords))
        selector_dims = map(dims, coords) do d, c
            _at_or_contains(d, c, atol)
        end
        selectors = map(val, DD.dims(selector_dims, DD.dims(x)))
        # Check the selector is in bounds / actually selectable
        I = DD.selectindices(DD.dims(x, dims), selectors; err=_False())::Union{Nothing,Tuple{Int,Int}}
        if !isnothing(I)
            D = map(rebuild, selector_dims, I)
            if x isa Raster
                prop = x[D]
                _ismissingval(missingval(x), prop) && return 0
                props = NamedTuple{K,Tuple{eltype(x)}}((prop,))
            else
                props = x[D][K]
                _ismissingval(missingval(x), props) && return 0
            end
            rows[i] = _maybe_add_fields(T, props, point, I)
            return 1
        end
    end
    return 0
end
# Normal point with missing / out of bounds data kept with `missing` fields
@inline function _extract_point!(
    row::Vector{T}, x::RasterStackOrArray, skipmissing::_False, point, i;
    dims=DD.dims(x, DEFAULT_POINT_ORDER),
    names::NamedTuple{K},
    atol=nothing,
    kw...
) where {T,K}
    # Get the actual dimensions available in the object
    coords = map(DD.commondims(x, dims)) do d
        _dimcoord(d, point)
    end
    # If any are coordinates missing, also return missing for everything
    rows[i] = if any(map(ismissing, coords))
        _maybe_add_fields(T, map(_ -> missing, names), missing, missing)
    else
        selector_dims = map(dims, coords) do d, c
            _at_or_contains(d, c, atol)
        end
        selectors = map(val, DD.dims(selector_dims, DD.dims(x)))
        # Check the selector is in bounds / actually selectable
        I = DD.selectindices(DD.dims(x, dims), selectors; err=_False())::Union{Nothing,Tuple{Int,Int}}
        if isnothing(I)
            _maybe_add_fields(T, map(_ -> missing, names), point, missing)
        else
            D = map(rebuild, selector_dims, I)
            props = if x isa Raster
                NamedTuple{K,Tuple{eltype(x)}}((x[D],))
            else
                NamedTuple(x[D])[K]
            end
            _maybe_add_fields(T, props, point, I)
        end
    end
    return 1
end

# Maybe add optional fields
# It is critically important for performance that this is type stable
Base.@assume_effects :total function _maybe_add_fields(
    ::Type{T}, props::NamedTuple, point, I
)::T where {T<:NamedTuple{K}} where K
    if :geometry in K
        :index in K ? merge((; geometry=point, index=I), props) : merge((; geometry=point), props)
    else
        :index in K ? merge((; index=I), props) : props
    end |> T
end
