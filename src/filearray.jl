
struct FileArray{X,T,N,F<:AbstractString,S<:Tuple,K,C} <: AbstractDiskArray{T,N}
    filename::F
    size::S
    key::K
    chunks::C
    write::Bool
end
function FileArray{X,T,N}(
    filename::F, size::S, key::K, chunks::C=nothing, write=false
) where {X,T,N,F,S<:Tuple,K,C}
    FileArray{X,T,N,F,S,K,C}(filename, size, key, chunks, write)
end
function FileArray{X,T,N}(filename::AbstractString, size::Tuple; 
    key=nothing, chunks=nothing, write=false
) where {X,T,N}
    FileArray{X,T,N}(filename, size, key, chunks, write)
end

filename(A::FileArray) = A.filename
key(A::FileArray) = A.key
chunks(A::FileArray) = A.chunks
Base.size(A::FileArray) = A.size

DiskArrays.haschunks(A::FileArray) = _haschunks(chunks(A))

_haschunks(::Nothing) = DiskArrays.Unchunked()
_haschunks(x) = DiskArrays.Chunked()

DiskArrays.readblock!(A::FileArray, dst, r::AbstractUnitRange...) = 
    open(o -> dst .= o[r...], A)
DiskArrays.writeblock!(A::FileArray, src, r::AbstractUnitRange...) = 
    open(o -> o[r...] .= src, A)

Base.open(f::Function, A::FileArray; kw...) = open(f, A; key=key(A), kw...)

function DiskArrays.eachchunk(f, A::FileArray) 
    open(A) do parent
        DiskArrays.GridChunks(parent, chunks(parent))
    end
end
