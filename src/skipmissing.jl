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

_missing(x, itr) = ismissing(x) || x == missingval(itr)  # mind the order, as comparison with missing returns missing
function _missing(x::AbstractFloat, itr)
    if isnothing(missingval(itr))
        return false
    elseif isnan(missingval(itr))
        return ismissing(x) || isnan(x)
    else
        return ismissing(x) || x == missingval(itr)
    end
end

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
