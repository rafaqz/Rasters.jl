"""
    OpenStack{X,K}

    OpenStack{X,K}(dataset)

A wrapper for any stack-like opened dataset that can be indexed
with `Symbol` keys to retrieve `AbstractArray` layers.

`OpenStack` is usually hidden from users, wrapped in a regular `RasterStack`
passed as the function argument in `open(stack)` when the stack is
contained in a single file.

`X` is a backend type like `NCDfile`, and `K` is a tuple of `Symbol` keys.
"""
struct OpenStack{X,K,DS}
    dataset::DS
end
OpenStack{X,K}(dataset::DS) where {X,K,DS} = OpenStack{X,K,DS}(dataset)

# OpenStack has `X` and `K` parameter that is not recoverable from fields.
ConstructionBase.constructorof(::Type{<:OpenStack{X,K}}) where {X,K} = OpenStack{X,K} 

dataset(os::OpenStack) = os.dataset
Base.keys(os::OpenStack{<:Any,K}) where K = K
# TODO test this, and does it make sense to return an iterator here?
Base.values(os::OpenStack{<:Any}) = (os[k] for k in keys(os))
# Indexing OpenStack returns memory-backed Raster. 
Base.getindex(os::OpenStack, key::Symbol) = dataset(os)[key]
