abstract type AbstractFileStack end

struct FileStack{X,K,F<:AbstractString,T<:NamedTuple,S<:NamedTuple,W,EC,HC} <: AbstractFileStack
    filename::F
    types::T
    sizes::S
    window::W
    eachchunk::EC
    haschunks::HC
    write::Bool
end
function FileStack{X,K}(
    filename::F, types::T, sizes::S, window::W, eachchunk::EC, haschunks::HC, write::Bool
) where {X,K,F,T,S,W,EC,HC}
    FileStack{X,K,F,T,S,W,EC,HC}(filename, types, sizes, window, eachchunk, haschunks, write)
end

ConstructionBase.constructorof(::Type{FileStack{X,K}}) where {X,K} = FileStack{X,K} 

filename(fs::FileStack) = fs.filename
window(fs::FileStack) = fs.window
Base.keys(fs::FileStack{<:Any,K}) where K = K
Base.values(fs::FileStack{<:Any}) = (fs[k] for k in keys(fs))
function Base.getindex(fs::FileStack{X}, key::Symbol) where X
    size = fs.sizes[key]
    eachchunk = fs.eachchunk[key]
    haschunks = fs.haschunks[key]
    T = fs.types[key]
    N = length(size)
    A = FileArray{X,T,N}(filename(fs), size, key, eachchunk, haschunks, fs.write)
    window(fs) isa Nothing ? A : view(A, window(fs)...)
end

function Base.view(st::FileStack, I...)
    @set st.window = reindex_window(st, st.window, I)
end
