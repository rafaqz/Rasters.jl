abstract type AbstractFileStack end

struct FileStack{X,K,F<:AbstractString,T<:NamedTuple,S<:NamedTuple} <: AbstractFileStack
    filename::F
    types::T
    sizes::S
    write::Bool
end
function FileStack{X,K}(filename::F, types::T, sizes::S, write::Bool) where {X,K,F,T,S}
    FileStack{X,K,F,T,S}(filename, types, sizes, write)
end

filename(fs::FileStack) = fs.filename
Base.keys(fs::FileStack{<:Any,K}) where K = K
Base.values(fs::FileStack{<:Any}) = (fs[k] for k in keys(fs))
function Base.getindex(fs::FileStack{X}, key::Symbol) where X
    size = fs.sizes[key]
    T = fs.types[key]
    N = length(size)
    FileArray{X,T,N}(filename(fs), size; key=key, write=fs.write)
end
