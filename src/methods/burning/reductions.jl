# _TakeFirst: `first` with a missing value check
# TODO: clarify why we need this  for `first` and not `last`
struct _TakeFirst{MV}
    missingval::MV
end
(tf::_TakeFirst)(a, b) = a === tf.missingval ? b : a
_take_last(a, b) = b

# Get the operator for a reduction. 
# This can be used for online stats style optimizations
_reduce_op(::typeof(sum)) = Base.add_sum
_reduce_op(::typeof(prod)) = Base.mul_prod
_reduce_op(::typeof(minimum)) = min
_reduce_op(::typeof(maximum)) = max
_reduce_op(::typeof(last)) = _take_last
_reduce_op(f, missingval) = _reduce_op(f)
_reduce_op(::typeof(first), missingval) = _TakeFirst(missingval)
_reduce_op(x) = nothing

# Can we split the reducion and combine later
_is_op_threadsafe(::typeof(sum)) = true
_is_op_threadsafe(::typeof(prod)) = true
_is_op_threadsafe(::typeof(minimum)) = true
_is_op_threadsafe(::typeof(maximum)) = true
_is_op_threadsafe(f) = false

# Get an initial value for the reduction
# TODO: clarify the difference between returning missing and zero
_reduce_init(reducer, st::AbstractRasterStack, missingval) = map(A -> _reduce_init(reducer, A, missingval), st)
_reduce_init(reducer, ::AbstractRaster{T}, missingval) where T = _reduce_init(reducer, T, missingval)
_reduce_init(reducer, nt::NamedTuple, missingval) = map(x -> _reduce_init(reducer, x, missingval), nt)
_reduce_init(f, x, missingval) = _reduce_init(f, typeof(x), missingval)
_reduce_init(::Nothing, x::Type{T}, missingval) where T = zero(T)
_reduce_init(f::Function, ::Type{T}, missingval) where T = zero(f((zero(nonmissingtype(T)), zero(nonmissingtype(T)))))
_reduce_init(::typeof(sum), ::Type{T}, missingval) where T = zero(nonmissingtype(T))
_reduce_init(::typeof(prod), ::Type{T}, missingval) where T = oneunit(nonmissingtype(T))
_reduce_init(::typeof(minimum), ::Type{T}, missingval) where T = typemax(nonmissingtype(T))
_reduce_init(::typeof(maximum), ::Type{T}, missingval) where T = typemin(nonmissingtype(T))
_reduce_init(::typeof(last), ::Type{T}, missingval) where T = _maybe_nothing_to_missing(missingval)