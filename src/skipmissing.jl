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


## This is adapted from https://github.com/JuliaLang/julia/blob/50713ee4a82eb1b5613647cd74b027315f665080/base/missing.jl#L273-L362
## It ensures more precise calculation of reductions like sum. See https://github.com/rafaqz/Rasters.jl/issues/812
Base.mapreduce(f, op, itr::SkipMissingVal) =
   Base._mapreduce(f, op, IndexStyle(itr.x), itr)

function Base._mapreduce(f, op, ::IndexLinear, itr::SkipMissingVal)
    A = itr.x
    ai = missing
    inds = LinearIndices(A)
    i = first(inds)
    ilast = last(inds)
    for outer i in i:ilast
        @inbounds ai = A[i]
        !_missing(ai, itr) && break
    end
    _missing(ai, itr) && return mapreduce_empty(f, op, eltype(itr))
    a1::eltype(itr) = ai
    i == typemax(typeof(i)) && return mapreduce_first(f, op, a1)
    i += 1
    ai = missing
    for outer i in i:ilast
        @inbounds ai = A[i]
        !_missing(ai, itr) && break
    end
    _missing(ai, itr) && return mapreduce_first(f, op, a1)
    # We know A contains at least two non-missing entries: the result cannot be nothing
    something(Base.mapreduce_impl(f, op, itr, first(inds), last(inds)))
end

Base._mapreduce(f, op, ::IndexCartesian, itr::SkipMissingVal) = Base.mapfoldl(f, op, itr)

Base.mapreduce_impl(f, op, A::SkipMissingVal, ifirst::Integer, ilast::Integer) =
    Base.mapreduce_impl(f, op, A, ifirst, ilast, Base.pairwise_blocksize(f, op))

# Returns nothing when the input contains only missing values, and Some(x) otherwise
@noinline function Base.mapreduce_impl(f, op, itr::SkipMissingVal,
                                  ifirst::Integer, ilast::Integer, blksize::Int)
    A = itr.x
    if ifirst > ilast
        return nothing
    elseif ifirst == ilast
        @inbounds a1 = A[ifirst]
        if ismissing(a1)
            return nothing
        else
            return Some(mapreduce_first(f, op, a1))
        end
    elseif ilast - ifirst < blksize
        # sequential portion
        ai = missing
        i = ifirst
        for outer i in i:ilast
            @inbounds ai = A[i]
            !_missing(ai, itr) && break
        end
        _missing(ai, itr) && return nothing
        a1 = ai::eltype(itr)
        i == typemax(typeof(i)) && return Some(mapreduce_first(f, op, a1))
        i += 1
        ai = missing
        for outer i in i:ilast
            @inbounds ai = A[i]
            !_missing(ai, itr) && break
        end
        _missing(ai, itr) && return Some(mapreduce_first(f, op, a1))
        a2 = ai::eltype(itr)
        i == typemax(typeof(i)) && return Some(op(f(a1), f(a2)))
        i += 1
        v = op(f(a1), f(a2))
        @simd for i = i:ilast
            @inbounds ai = A[i]
            if !_missing(ai, itr)
                v = op(v, f(ai))
            end
        end
        return Some(v)
    else
        # pairwise portion
        imid = ifirst + (ilast - ifirst) >> 1
        v1 = Base.mapreduce_impl(f, op, itr, ifirst, imid, blksize)
        v2 = Base.mapreduce_impl(f, op, itr, imid+1, ilast, blksize)
        if v1 === nothing && v2 === nothing
            return nothing
        elseif v1 === nothing
            return v2
        elseif v2 === nothing
            return v1
        else
            return Some(op(something(v1), something(v2)))
        end
    end
end