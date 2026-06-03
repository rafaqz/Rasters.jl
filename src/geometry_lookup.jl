"""
    GeometryLookup(data, dims = (X(), Y()); geometrycolumn = nothing)

A lookup type for geometry dimensions in vector data cubes.

`GeometryLookup` provides efficient spatial indexing and lookup for
geometries using an STRtree (Sort-Tile-Recursive tree).

It is used as the lookup type for geometry dimensions in vector
data cubes, enabling fast spatial queries and operations.

It spans the dimensions given to it in `dims`, as well as the dimension
 it's wrapped in - you would construct a DimArray with a GeometryLookup
like `DimArray(data, Geometry(GeometryLookup(data, dims)))`.
Here, `Geometry` is a dimension - but selectors in X and Y will also work!

# Examples

```julia
using Rasters

using NaturalEarth
import GeometryOps as GO

# construct the polygon lookup
polygons = NaturalEarth.naturalearth("admin_0_countries", 110).geometry
polygon_lookup = GeometryLookup(polygons, (X(), Y()))

# create a DimArray with the polygon lookup
dv = rand(Geometry(polygon_lookup))

# select the polygon with the centroid of the 88th polygon
dv[Geometry(Contains(GO.centroid(polygons[88])))] == dv[Geometry(88)] # true
```
"""
struct GeometryLookup{T,A<:AbstractVector{T},D,Tree,CRS,MD} <: DD.Dimensions.MultiDimensionalLookup{T}
    data::A
    tree::Tree
    dims::D
    crs::CRS
    metadata::MD
end
function GeometryLookup(data, dims=(X(), Y());
    geometrycolumn=nothing, crs=nokw, tree=nokw, metadata=NoMetadata(),
)
    # First, retrieve the geometries - from a table, vector of geometries, etc.
    geometries = _get_geometries(data, geometrycolumn)
    geometries = Missings.disallowmissing(geometries)

    if isnokw(crs)
        crs = GI.crs(data)
        if isnothing(crs)
            crs = GI.crs(first(geometries))
        end
    end

    if any(!GI.isgeometry, geometries)
        throw(ArgumentError("""
        The collection passed in to `GeometryLookup` has some elements that are not geometries
        (`GI.isgeometry(x) == false` for some `x` in `data`).
        """))
    end
    if length(dims) != 2
        throw(ArgumentError("""
        The `dims` argument to `GeometryLookup` must have two dimensions, but it has $(length(dims)) dimensions (`$(dims)`).
        Please make sure that it has only two dimensions, like `(X(), Y())`.
        """))
    end
    tree = _resolve_geometry_tree(tree, geometries)
    GeometryLookup(geometries, tree, dims, crs, metadata)
end

# Build/accept the spatial index tree for `GeometryLookup`.
function _resolve_geometry_tree(tree, geometries)
    isnokw(tree) && return SortTileRecursiveTree.STRtree(geometries)
    isnothing(tree) && return nothing
    if GO.SpatialTreeInterface.isspatialtree(tree)
        return tree isa DataType ? tree(geometries) : tree
    end
    throw(ArgumentError("""
    Got an argument for `tree` which is not a valid spatial tree (according to `GeometryOps.SpatialTreeInterface`)
    nor `nokw` or `nothing`

    Type is $(typeof(tree))
    """))
end

GeoInterface.crs(l::GeometryLookup) = l.crs
setcrs(l::GeometryLookup, crs) = rebuild(l; crs)

# DimensionalData protocol --------------------------------------------------

DD.dims(l::GeometryLookup) = l.dims
# Expose the inner X/Y dims so DD's dim resolution can descend through the
# Geometry dim - matches what DD does for MergedLookup (see merged.jl:70).
# Without this, `_dim_query1(f, op, t, geom_dim, query)` falls through to the
# `dims(x)` branch in primitives.jl:796 which returns the dim itself,
# causing infinite recursion on any partial query (e.g. `a[Ti=1]` when the
# array also has a Geometry dim).
DD.dims(d::DD.Dimension{<:GeometryLookup}) = dims(val(d))
DD.order(::GeometryLookup) = Lookups.Unordered()
DD.parent(lookup::GeometryLookup) = lookup.data
DD.metadata(l::GeometryLookup) = l.metadata
# TODO: format for geometry lookup
DD.Dimensions.format(l::GeometryLookup, D::Type, values, axis::AbstractRange) = l

# Rebuild the tree when data changes, unless the caller passes one in.
function DD.rebuild(lookup::GeometryLookup;
    data=lookup.data, tree=nokw,
    dims=lookup.dims, crs=nokw,
    metadata=lookup.metadata,
)
    new_tree = if isnokw(tree)
        if data == lookup.data
            lookup.tree
        elseif isempty(data)
            nothing
        else
            SortTileRecursiveTree.STRtree(data)
        end
    else
        _resolve_geometry_tree(tree, data)
    end
    new_crs = if isnokw(crs)
        data_crs = GI.crs(data)
        isnothing(data_crs) ? lookup.crs : data_crs
    else
        crs
    end
    return GeometryLookup(Missings.disallowmissing(data), new_tree, dims, new_crs, metadata)
end

function Lookups.bounds(lookup::GeometryLookup)
    if isempty(lookup.data)
        return Extents.Extent(NamedTuple{DD.name.(lookup.dims)}(ntuple(_ -> (nothing, nothing), 2)))
    end
    if isnothing(lookup.tree)
        return mapreduce(GI.extent, Extents.union, lookup.data)
    end
    return GI.extent(lookup.tree)
end

# Return an `Int` or `Vector{Bool}` for index selection.
# Base: a standard index can go straight into getindex.
Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.StandardIndices) = sel
function Lookups.selectindices(lookup::GeometryLookup, sel::DD.DimTuple)
    selectindices(lookup, map(_val_or_nothing, sortdims(sel, dims(lookup))))
end
function Lookups.selectindices(lookup::GeometryLookup, sel::NamedTuple{K}) where K
    dimsel = map(rebuild, map(name2dim, K), values(sel))
    selectindices(lookup, dimsel)
end
function Lookups.selectindices(lookup::GeometryLookup, sel::Tuple)
    if (length(sel) == length(dims(lookup))) && all(map(s -> s isa At, sel))
        i = findfirst(x -> all(map(Dimensions._matches, sel, x)), lookup)
        isnothing(i) && _coord_not_found_error(sel)
        return i
    else
        return [Dimensions._matches(sel, x) for x in lookup]
    end
end
function Lookups.selectindices(lookup::GeometryLookup, sel::Contains)
    sel_ext = GI.extent(val(sel))
    potential_candidates = _maybe_get_candidates(lookup, sel_ext)
    filter(potential_candidates) do candidate
        GO.contains(lookup.data[candidate], val(sel))
    end
end
function Lookups.selectindices(lookup::GeometryLookup, sel::At)
    geom = val(sel)
    @assert GI.isgeometry(geom)
    candidates = _maybe_get_candidates(lookup, GI.extent(geom))
    x = findfirst(candidates) do candidate
        GO.equals(geom, lookup.data[candidate])
    end
    isnothing(x) && throw(ArgumentError("$sel not found in lookup"))
    return x
end
function Lookups.selectindices(lookup::GeometryLookup, sel::Near)
    geom = val(sel)
    @assert GI.isgeometry(geom)
    # TODO: temporary
    @assert GI.trait(geom) isa GI.PointTrait "Only point geometries are supported for the near lookup at this point!  We will add more geometry support in the future."
    # TODO: branch and bound on the spatial tree instead of a full scan.
    return findmin(x -> GO.distance(geom, x), lookup.data)[2]
end
function Lookups.selectindices(lookup::GeometryLookup, sel::Touches)
    sel_ext = GI.extent(val(sel))
    potential_candidates = _maybe_get_candidates(lookup, sel_ext)
    return filter(potential_candidates) do candidate
        GO.intersects(lookup.data[candidate], val(sel))
    end
end
function Lookups.selectindices(lookup::GeometryLookup,
    (xs, ys)::Tuple{Union{<:Touches},Union{<:Touches}},
)
    target_ext = Extents.Extent(X=(first(xs), last(xs)), Y=(first(ys), last(ys)))
    potential_candidates = _maybe_get_candidates(lookup, target_ext)
    return filter(potential_candidates) do candidate
        GO.intersects(lookup.data[candidate], target_ext)
    end
end
function Lookups.selectindices(lookup::GeometryLookup,
    (xs, ys)::Tuple{Union{<:DD.IntervalSets.ClosedInterval},Union{<:DD.IntervalSets.ClosedInterval}},
)
    target_ext = Extents.Extent(X=extrema(xs), Y=extrema(ys))
    potential_candidates = _maybe_get_candidates(lookup, target_ext)
    filter(potential_candidates) do candidate
        GO.covers(target_ext, lookup.data[candidate])
    end
end
function Lookups.selectindices(lookup::GeometryLookup,
    (x, y)::Tuple{Union{<:At,<:Contains},Union{<:At,<:Contains}},
)
    xval, yval = val(x), val(y)
    lookup_ext = Lookups.bounds(lookup)
    if lookup_ext.X[1] <= xval <= lookup_ext.X[2] && lookup_ext.Y[1] <= yval <= lookup_ext.Y[2]
        potential_candidates = GO.SpatialTreeInterface.query(lookup.tree, (xval, yval))
        isempty(potential_candidates) && return Int[]
        filter(potential_candidates) do candidate
            GO.contains(lookup.data[candidate], (xval, yval))
        end
    else
        return Int[]
    end
end

@inline Lookups.reducelookup(l::GeometryLookup) = NoLookup(OneTo(1))

function Lookups.show_properties(io::IO, mime, lookup::GeometryLookup)
    print(io, " ")
    show(IOContext(io, :inset => "", :dimcolor => 244), mime, DD.basedims(lookup))
end

# Dimension methods ---------------------------------------------------------

@inline _reducedims(lookup::GeometryLookup, dim::DD.Dimension) =
    rebuild(dim, [map(x -> zero(x), dim.val[1])])

function DD.format(dim::DD.Dimension{<:GeometryLookup}, axis::AbstractRange)
    checkaxis(dim, axis)
    return dim
end

# Helpers -------------------------------------------------------------------

_val_or_nothing(::Nothing) = nothing
_val_or_nothing(d::DD.Dimension) = val(d)

# Get the candidates for the selector extent. If the selector extent is
# disjoint from the tree rootnode extent, return an empty index.
function _maybe_get_candidates(lookup::GeometryLookup, selector_extent)
    tree = lookup.tree
    isnothing(tree) && return 1:length(lookup)
    Extents.disjoint(GI.extent(tree), selector_extent) && return Int[]
    potential_candidates = GO.SpatialTreeInterface.query(
        tree,
        Base.Fix1(Extents.intersects, selector_extent)
    )
    isempty(potential_candidates) && return Int[]
    return potential_candidates
end
