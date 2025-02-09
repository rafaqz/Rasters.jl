"""
    OpenStack{X,K}

    OpenStack{X,K}(dataset)

A wrapper for any stack-like opened dataset that can be indexed
with `Symbol` keys to retrieve `AbstractArray` layers.

`OpenStack` is usually hidden from users, wrapped in a regular `RasterStack`
passed as the function argument in `open(stack)` when the stack is
contained in a single file.

`X` is a backend type like `NCDsource`, and `K` is a tuple of `Symbol` keys.
"""
struct OpenStack{X,K,T,DS,M}
    dataset::DS
    mods::M
end
function OpenStack{X,K,T}(
    dataset::DS, mods::M=NoMod()
) where {X,K,T,DS,M}
    OpenStack{X,K,T,DS,M}(dataset, mods)
end

dataset(os::OpenStack) = os.dataset

# OpenStack has `X` and `K` parameter that is not recoverable from fields.
ConstructionBase.constructorof(::Type{<:OpenStack{X,K,T}}) where {X,K,T} = 
    OpenStack{X,K,T} 

DD.data_eltype(::OpenStack{<:Any,<:Any,T}) where T = T

Base.keys(::OpenStack{<:Any,K}) where K = K
# TODO test this, and does it make sense to return an iterator here?
Base.values(os::OpenStack{<:Any,K}) where K = (os[k] for k in K)
# Indexing OpenStack returns memory-backed Raster. 
function Base.getindex(os::OpenStack{<:Any,K}, key::Symbol) where K
    mods = os.mods
    if mods isa AbstractModifications
        _maybe_modify(dataset(os)[key], mods)
    else
        _maybe_modify(dataset(os)[key], NamedTuple{K}(mods)[key])
    end
end
