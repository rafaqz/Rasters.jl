
"""
    FileArray{X} <: DiskArrays.AbstractDiskArray

Filearray is a DiskArrays.jl `AbstractDiskArray`. Instead of holding
an open object, it just holds a filename string that is opened lazily 
when it needs to be read.
"""
struct FileArray{X,T,N,K,EC,HC} <: DiskArrays.AbstractDiskArray{T,N}
    filename::String
    size::NTuple{N,Int}
    key::K
    eachchunk::EC
    haschunks::HC
    write::Bool
end
function FileArray{X,T,N}(
    filename, size, key::K, eachchunk::EC=size, 
    haschunks::HC=DA.Unchunked(), write=false
) where {X,T,N,K,EC,HC}
    FileArray{X,T,N,K,EC,HC}(filename, size, key, eachchunk, haschunks, write)
end
function FileArray{X,T,N}(filename::String, size::Tuple; 
    key=nothing, eachchunk=size, haschunks=DA.Unchunked(), write=false
) where {X,T,N}
    FileArray{X,T,N}(filename, size, key, eachchunk, haschunks, write)
end

# FileArray has X, T and N parameters not recoverable from fields
ConstructionBase.constructorof(::Type{<:FileArray{X,T,N}}) where {X,T,N} = FileArray{X,T,N}

filename(A::FileArray) = A.filename
key(A::FileArray) = A.key
Base.size(A::FileArray) = A.size
DA.eachchunk(A::FileArray) = A.eachchunk
DA.haschunks(A::FileArray) = A.haschunks

# Run function `f` on the result of _open for the file type
function Base.open(f::Function, A::FileArray{X}; write=A.write, kw...) where X
    _open(f, X, filename(A); key=key(A), write, kw...)
end

function DA.readblock!(A::FileArray, dst, r::AbstractUnitRange...)
    open(A) do O
        DA.readblock!(O, dst, r...)
    end
end
function DA.writeblock!(A::FileArray, src, r::AbstractUnitRange...) 
    open(A; write=A.write) do O
        DA.writeblock!(O, src, r...)
    end
end


"""
    RasterDiskArray <: DiskArrays.AbstractDiskArray

A basic DiskArrays.jl wrapper for objects that don't have one defined yet. 
When we `open` a `FileArray` it is replaced with a `RasterDiskArray`.
"""
struct RasterDiskArray{X,T,N,V,EC,HC} <: DiskArrays.AbstractDiskArray{T,N}
    var::V
    eachchunk::EC
    haschunks::HC
end
function RasterDiskArray{X}(
    var::V, eachchunk=DA.eachchunk(var), haschunks=DA.haschunks(var)
) where {X,V}
    T = eltype(var)
    N = ndims(var)
    RasterDiskArray{X,T,N,V,typeof(eachchunk),typeof(haschunks)}(var, eachchunk, haschunks)
end

Base.parent(A::RasterDiskArray) = A.var
Base.size(A::RasterDiskArray{<:Any,T,N}) where {T,N} = size(parent(A))::NTuple{N,Int}

DA.haschunks(A::RasterDiskArray) = A.haschunks
DA.eachchunk(A::RasterDiskArray) = A.eachchunk
DA.readblock!(A::RasterDiskArray, aout, r::AbstractUnitRange...) = aout .= parent(A)[r...]
DA.writeblock!(A::RasterDiskArray, v, r::AbstractUnitRange...) = parent(A)[r...] .= v

struct MissingDiskArray{T,N,V} <: DiskArrays.AbstractDiskArray{T,N}
    var::V
end
function MissingDiskArray(::Type{MT}, var::A) where {MT,A <: AbstractArray{T,N}} where {T,N}
    MissingDiskArray{MT,N,A}(var)
end

Base.parent(A::MissingDiskArray) = A.var
Base.size(A::MissingDiskArray) = size(parent(A))

DA.haschunks(A::MissingDiskArray) = DA.haschunks(parent(A))
DA.eachchunk(A::MissingDiskArray) = DA.eachchunk(parent(A))
DA.readblock!(A::MissingDiskArray, aout, r::AbstractUnitRange...) = DA.readblock!(parent(A), aout, r...)
DA.writeblock!(A::MissingDiskArray, v, r::AbstractUnitRange...) = DA.writeblock!(parent(A), v, r...)
