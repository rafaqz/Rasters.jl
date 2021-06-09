abstract type AbstractFileStack end

struct FileStack{X,K,F<:AbstractString,D<:Tuple,T<:NamedTuple,S<:NamedTuple,W,EC,HC} <: AbstractFileStack
    filename::F
    dims::D
    types::T
    sizes::S
    window::W
    eachchunk::EC
    haschunks::HC
    write::Bool
end
function FileStack{X,K}(
    filename::F, dims::D, types::T, sizes::S, window::W, eachchunk::EC, haschunks::HC, write::Bool
) where {X,K,F,D,T,S,W,EC,HC}
    FileStack{X,K,F,D,T,S,W,EC,HC}(filename, dims, types, sizes, window, eachchunk, haschunks, write)
end

ConstructionBase.constructorof(::Type{<:FileStack{X,K}}) where {X,K} = FileStack{X,K} 

filename(fs::FileStack) = fs.filename
window(fs::FileStack) = fs.window
dims(fs::FileStack) = fs.dims
window(fs::NamedTuple) = nothing
Base.keys(fs::FileStack{<:Any,K}) where K = K
Base.values(fs::FileStack{<:Any}) = (fs[k] for k in keys(fs))
function Base.getindex(fs::FileStack{X}, key::Symbol) where X
    size = fs.sizes[key]
    eachchunk = fs.eachchunk[key]
    haschunks = fs.haschunks[key]
    T = fs.types[key]
    N = length(size)
    A = FileArray{X,T,N}(filename(fs), size, key, eachchunk, haschunks, fs.write)
end

function Base.view(st::FileStack, I...)
end

@propagate_inbounds function Base.view(fs::FileStack, I...)
    I = dims2indices(fs, I...)
    window, refwindow = slicedims(view, window(fs), I)
    @set filestack.window = window
    rebuild(s, filestack)
end
