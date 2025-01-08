"""
    GeometryLookup(data, dims = (X(), Y()); geometrycolumn = nothing)


The other thing I'm thinking of is that we could have a feature collection in `data` so that you can persist per geom attrs

A lookup type for geometry dimensions in vector data cubes.

`GeometryLookup` provides efficient spatial indexing and lookup for geometries using an STRtree (Sort-Tile-Recursive tree).
It is used as the lookup type for geometry dimensions in vector data cubes, enabling fast spatial queries and operations.

It spans the dimensions given to it in `dims`, as well as the dimension it's wrapped in - you would construct a DimArray with a GeometryLookup
like `DimArray(data, Geometry(GeometryLookup(data, dims)))`.  Here, `Geometry` is a dimension - but selectors in X and Y will also eventually work!


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
struct GeometryLookup{T, D} <: Lookups.Lookup{T, 1}
    data::Vector{T}
    tree::SortTileRecursiveTree.STRtree
    dims::D
end

function GeometryLookup(data, dims = (X(), Y()); geometrycolumn = nothing)
    # First, retrieve the geometries - from a table, vector of geometries, etc.
    geometries = _get_geometries(data, geometrycolumn)
    # Check that the geometries are actually geometries
    if any(!GI.isgeometry, geometries)
        throw(ArgumentError("""
        The collection passed in to `GeometryLookup` has some elements that are not geometries 
        (`GI.isgeometry(x) == false` for some `x` in `data`).
        """))
    end
    # Make sure there are only two dimensions
    if length(dims) != 2
        throw(ArgumentError("""
        The `dims` argument to `GeometryLookup` must have two dimensions, but it has $(length(dims)) dimensions (`$(dims)`).
        Please make sure that it has only two dimensions, like `(X(), Y())`.
        """))
    end
    # Build the lookup accelerator tree
    tree = SortTileRecursiveTree.STRtree(geometries)
    GeometryLookup(geometries, tree, dims)
end

DD.dims(l::GeometryLookup) = l.dims
DD.dims(d::DD.Dimension{<: GeometryLookup}) = val(d).dims
 
DD.order(::GeometryLookup) = Lookups.Unordered()

DD.parent(lookup::GeometryLookup) = lookup.data

DD.Dimensions.format(l::GeometryLookup, D::Type, values, axis::AbstractRange) = l

function DD.rebuild(lookup::GeometryLookup; data = lookup.data, tree = nothing, dims = lookup.dims)
    new_tree = if data == lookup.data
        lookup.tree
    else
        SortTileRecursiveTree.STRtree(data)
    end
    GeometryLookup(data, new_tree, dims)
end

# Return an `Int` or  Vector{Bool}
Lookups.selectindices(lookup::GeometryLookup, sel::DD.DimTuple) =
    selectindices(lookup, map(_val_or_nothing, sortdims(sel, dims(lookup))))
function Lookups.selectindices(lookup::GeometryLookup, sel::NamedTuple{K}) where K
    dimsel = map(rebuild, map(name2dim, K), values(sel))
    selectindices(lookup, dimsel) 
end
Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.StandardIndices) = sel
function Lookups.selectindices(lookup::GeometryLookup, sel::Tuple)
    if (length(sel) == length(dims(lookup))) && all(map(s -> s isa At, sel))
        i = findfirst(x -> all(map(_matches, sel, x)), lookup)
        isnothing(i) && _coord_not_found_error(sel)
        return i
    else
        return [_matches(sel, x) for x in lookup]
    end
end

function _maybe_get_candidates(tree, selector_extent)
    if Extents.disjoint(tree.rootnode.extent, selector_extent)
        return error("""
            The geometry with extent $(GI.extent(val(sel))) is outside of the extent of the geometry lookup.

            The geometry lookup has extent $(GI.extent(tree.rootnode.extent))
        """)
    end
    potential_candidates = SortTileRecursiveTree.query(tree, selector_extent)

    isempty(potential_candidates) && return error("""
        The geometry with extent $(GI.extent(val(sel))) does not interact with any of the geometries in the lookup.
    """)

    return potential_candidates
end

function Lookups.selectindices(lookup::GeometryLookup, sel::Contains)
    lookup_ext = lookup.tree.rootnode.extent
    sel_ext = GI.extent(val(sel))
    potential_candidates = _maybe_get_candidates(lookup.tree, sel_ext)

    for candidate in potential_candidates
        if GO.contains(lookup.data[candidate], val(sel))
            return candidate
        end
    end
    return error("""
        The geometry with extent $(GI.extent(val(sel))) is not contained by any of the geometries in the lookup.
    """)
end


function Lookups.selectindices(lookup::GeometryLookup, sel::Touches)
    lookup_ext = lookup.tree.rootnode.extent
    sel_ext = GI.extent(val(sel))
    potential_candidates = _maybe_get_candidates(lookup.tree, sel_ext)

    for candidate in potential_candidates
        if GO.intersects(lookup.data[candidate], val(sel))
            return candidate
        end
    end
    return error("""
        The geometry with extent $(GI.extent(val(sel))) does not touch any of the geometries in the lookup.
    """)
end

function Lookups.selectindices(lookup::GeometryLookup, (xs, ys)::Tuple{Union{<: At, <: Contains}, Union{<: At, <: Contains}})
    xval, yval = val(xs), val(ys)

    lookup_ext = lookup.tree.rootnode.extent

    if lookup_ext.X[1] <= xval <= lookup_ext.X[2] && lookup_ext.Y[1] <= yval <= lookup_ext.Y[2]
        potential_candidates = SortTileRecursiveTree.query(lookup.tree, (xval, yval))
        if isempty(potential_candidates)
            return error("""
            The point ($xval, $yval) is not within any of the geometries in the lookup.
        """) # no geometry intersects with it
        else
            for candidate in potential_candidates
                if GO.contains(lookup.data[candidate], (xval, yval))
                    return candidate
                end
            end
            return error("""
                The point ($xval, $yval) is not within any of the geometries in the lookup.
            """) # no geometry intersects with it
        end
    else
        return  error("""
        The point ($xval, $yval) is outside of the extent of the geometry lookup. 
        """) # outside of extent / geometrylookup bbox
    end
end

@inline Lookups.reducelookup(l::GeometryLookup) = NoLookup(OneTo(1))

function Lookups.show_properties(io::IO, mime, lookup::GeometryLookup)
    print(io, " ")
    show(IOContext(io, :inset => "", :dimcolor => 244), mime, DD.basedims(lookup))
end

# Dimension methods

@inline _reducedims(lookup::GeometryLookup, dim::DD.Dimension) =
    rebuild(dim, [map(x -> zero(x), dim.val[1])])

function DD._format(dim::DD.Dimension{<:GeometryLookup}, axis::AbstractRange)
    checkaxis(dim, axis)
    return dim
end

# Local functions
_val_or_nothing(::Nothing) = nothing
_val_or_nothing(d::DD.Dimension) = val(d)

#=



DD.@dim Geometry

using NaturalEarth
polygons = NaturalEarth.naturalearth("admin_0_countries", 110).geometry

polygon_lookup = GeometryLookup(polygons, SortTileRecursiveTree.STRtree(polygons), (X(), Y()))

dv = rand(Geometry(polygon_lookup))


@test_throws ErrorException dv[Geometry=(X(At(1)), Y(At(2)))]
@test dv[Geometry(Contains(GO.centroid(polygons[88])))] == dv[Geometry(88)]





Geometry(Where(GO.contains(geom2)))

X(At(1)), Y(At(2)), Ti(Near(2010))


# TODO: 
# metadata for crs
# 
=#