
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

_missing(x, itr) = x === missingval(itr) || x === missing

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
