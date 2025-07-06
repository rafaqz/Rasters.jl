#=
# Lookups methods

Here we define the methods for the Lookups API.
The main entry point is `selectindices`, which is used to select the indices of the geometries that match the selector.

We need to define methods that take selectors and convert them to extents, then GeometryOps needs 
=#

# Bounds - get the bounds of the lookup
Lookups.bounds(lookup::GeometryLookup) = if isempty(lookup.data)
    Extents.Extent(NamedTuple{DD.name.(lookup.dims)}(ntuple(2) do i; (nothing, nothing); end))
else
    if isnothing(lookup.tree)
        mapreduce(GI.extent, Extents.union, lookup.data)
    else
        GI.extent(lookup.tree)
    end
end

# Return an `Int` or Vector{Bool}
# Base case: got a standard index that can go into getindex on a base Array
Lookups.selectindices(lookup::GeometryLookup, sel::Lookups.StandardIndices) = sel
# other cases: 
# - decompose selectors
function Lookups.selectindices(lookup::GeometryLookup, sel::DD.DimTuple)
    selectindices(lookup, map(_val_or_nothing, sortdims(sel, dims(lookup))))
end
function Lookups.selectindices(lookup::GeometryLookup, sel::NamedTuple{K}) where K
    dimsel = map(rebuild, map(name2dim, K), values(sel))
    selectindices(lookup, dimsel) 
end
function Lookups.selectindices(lookup::GeometryLookup, sel::Tuple)
    if (length(sel) == length(dims(lookup))) && all(map(s -> s isa At, sel))
        i = findfirst(x -> all(map(Lookups._matches, sel, x)), lookup)
        isnothing(i) && _coord_not_found_error(sel)
        return i
    else
        return [Lookups._matches(sel, x) for x in lookup]
    end
end

# Selector implementations that use geometry (like Contains, Touches, etc.)
"""
    _maybe_get_candidates(lookup, selector_extent)

Get the candidates for the selector extent.  If the selector extent is disjoint from the tree rootnode extent, return an error.
"""
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

function Lookups.selectindices(lookup::GeometryLookup, sel::Contains)
    potential_candidates = _maybe_get_candidates(lookup, sel_ext)
    filter(potential_candidates) do candidate
        GO.contains(lookup.data[candidate], val(sel))
    end
end

function Lookups.selectindices(lookup::GeometryLookup, sel::At)
    if GI.trait(val(sel)) isa GI.PointTrait
        Lookups.selectindices(lookup, (At(GI.x(val(sel))), At(GI.y(val(sel)))))
    else # invoke the default method
        Lookups.at(lookup, sel)
    end
end

function Lookups.selectindices(lookup::GeometryLookup, sel::Near)
    geom = val(sel)
    @assert GI.isgeometry(geom)
    # TODO: temporary
    @assert GI.trait(geom) isa GI.PointTrait "Only point geometries are supported for the near lookup at this point!  We will add more geometry support in the future."

    # Get the nearest geometry
    # TODO: this sucks!  Use some branch and bound algorithm
    # on the spatial tree instead.
    # if pointtrait
    return findmin(x -> GO.distance(geom, x), lookup.data)[2]
    # else
    #     findmin(x -> GO.distance(GO.GEOS(), geom, x), lookup.data)[2]
    # end 
    # this depends on LibGEOS being installed.

end


function Lookups.selectindices(lookup::GeometryLookup, sel::Touches)
    sel_ext = GI.extent(val(sel))
    potential_candidates = _maybe_get_candidates(lookup, sel_ext)
    return filter(potential_candidates) do candidate
        GO.intersects(lookup.data[candidate], val(sel))
    end
end

function Lookups.selectindices(lookup::GeometryLookup, (xs, ys)::Tuple{Union{ <: Touches}, Union{ <: Touches}})
    target_ext = Extents.Extent(X = (first(xs), last(xs)), Y = (first(ys), last(ys)))
    potential_candidates = _maybe_get_candidates(lookup, target_ext)
    return filter(potential_candidates) do candidate
        GO.intersects(lookup.data[candidate], target_ext)
    end
end


function Lookups.selectindices(lookup::GeometryLookup, (xs, ys)::Tuple{Union{<: DD.IntervalSets.ClosedInterval}, Union{<: DD.IntervalSets.ClosedInterval}})
    target_ext = Extents.Extent(X = extrema(xs), Y = extrema(ys))
    potential_candidates = _maybe_get_candidates(lookup, target_ext)
    filter(potential_candidates) do candidate
        GO.covers(target_ext, lookup.data[candidate])
    end
end

function Lookups.selectindices(lookup::GeometryLookup, (xs, ys)::Tuple{Union{<: At, <: Contains, <: Real}, Union{<: At, <: Contains, <: Real}})
    xval, yval = val(xs), val(ys)

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


