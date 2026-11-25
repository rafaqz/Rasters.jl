"""
    skipmissing(itr::Raster)

Returns an iterable over the elements in a `Raster` object, skipping any values equal to either the `missingval` or `missing`.
"""
function Base.skipmissing(itr::Raster)
    if ismissing(missingval(itr))
        Base.SkipMissing(itr)
    else
        SkipMissingVal(itr)
    end
end

struct SkipMissingVal{T}
    x::T
end
Base.IteratorSize(::Type{<:SkipMissingVal}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{SkipMissingVal{T}}) where {T} = Base.IteratorEltype(T)
Base.eltype(::Type{SkipMissingVal{T}}) where {T} = Base.nonmissingtype(eltype(T))
missingval(itr::SkipMissingVal) = missingval(itr.x)

function Base.iterate(itr::SkipMissingVal, state...)
    y = iterate(itr.x, state...)
    y === nothing && return nothing
    item, state = y
    # We check for both `missing` and the raster `missingval`
    # Mostly the compiler should elide the `missing` check?
    while _missing(item, itr)
        y = iterate(itr.x, state)
        y === nothing && return nothing
        item, state = y
    end
    item, state
end

_missing(x, itr) = isequal(x, missingval(itr))
_missing(x::Missing, itr) = true
_missing(x::Nothing, itr) = false

Base.IndexStyle(::Type{<:SkipMissingVal{T}}) where {T} = IndexStyle(T)
Base.eachindex(itr::SkipMissingVal) =
    Iterators.filter(i -> !_missing(@inbounds(itr.x[i]), itr), eachindex(itr.x))
Base.keys(itr::SkipMissingVal) =
    Iterators.filter(i -> !_missing(@inbounds(itr.x[i]), itr), keys(itr.x))
@propagate_inbounds function Base.getindex(itr::SkipMissingVal, I...)
    v = itr.x[I...]
    _missing(v, itr) && throw(MissingException("the value at index $I is the raster missingval"))
    v
end

function Base.show(io::IO, s::SkipMissingVal)
    print(io, "skipmissing(")
    show(io, s.x)
    print(io, ')')
end


# Specialize broadcasted passmissing to respect raster.
struct _RasterPassMissing{F, M, O} <: Function
    f::F
    missingvals::M  # missingvals of each broadcast arg
    out_mv::O       # missingval of resulting raster
end

@inline function (pm::_RasterPassMissing)(xs...)
    for (x, mv) in zip(xs, pm.missingvals)
        ismissing(x) && return pm.out_mv
        !isnothing(mv) && x === mv && return pm.out_mv
    end
    return pm.f(xs...)
end

# Specialize broadcasted calls with a function wrapped in PassMissing
function Base.Broadcast.broadcasted(f::Missings.PassMissing, A::AbstractRaster, args...)
    all_args = (A, args...)
    mvs = missingval.(all_args)
    mv = missingval(A)
    # Attempt to fix the missingval before any broadcast
    ElType = Base.Broadcast.combine_eltypes(f.f, all_args)
    if Base.isconcretetype(ElType) # without a concrete type this won't work
        mv_fixed = _fix_missingval(ElType, mv)
        mv = isnothing(mv_fixed) ? mv : mv_fixed
    end
    new_f = _RasterPassMissing(f.f, mvs, mv)
    return Base.Broadcast.broadcasted(new_f, A, args...)
end