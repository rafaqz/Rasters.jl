struct LazySubArray{T,I}
    data::T
    dims::I
end

Base.parent(a::LazySubArray) = a.data

dims(a::LazySubArray) = a.dims

# Rebuild the stack without the LazySubArray, and index with the LazySubArray indices
@inline Base.getindex(stack::AbstractGeoStack{LazySubArray}, key::Key) = begin
    lsa = parent(stack)
    data(rebuild(stack, parent(lsa)), key, dims2indices(dims(lsa))...)
end
