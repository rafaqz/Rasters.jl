"""
    FileArray{S} <: DiskArrays.AbstractDiskArray

Filearray is a DiskArrays.jl `AbstractDiskArray`. Instead of holding
an open object, it just holds a filename string that is opened lazily
when it needs to be read.
"""
struct FileArray{S,T,N,Na,G,EC,HC,M<:AbstractModifications} <: DiskArrays.AbstractDiskArray{T,N}
    filename::String
    size::NTuple{N,Int}
    name::Na
    group::G
    eachchunk::EC
    haschunks::HC
    mod::M
    write::Bool
end
function FileArray{S,T,N}(
    filename,
    size::NTuple{N,Int},
    name::Na,
    group::G,
    eachchunk::EC,
    haschunks::HC,
    mod::M,
    write::Bool,
) where {S,T,N,Na,G,EC,HC,M}
    FileArray{S,T,N,Na,G,EC,HC,M}(
        String(filename), size, name, group, eachchunk, haschunks, mod, write
    )
end
function FileArray{S,T,N}(filename::AbstractString, size::Tuple;
    name=nokw,
    group=nokw,
    eachchunk=size,
    haschunks=DA.Unchunked(),
    mod,
    write=false
) where {S,T,N}
    name = isnokw(name) ? nothing : name
    group = isnokw(group) ? nothing : group
    FileArray{S,T,N}(filename, size, name, group, eachchunk, haschunks, mod, write)
end
function FileArray{S}(
    var::AbstractArray{<:Any,N}, filename; mod, kw...
) where {S,N}
    eachchunk = DA.eachchunk(var)
    haschunks = DA.haschunks(var)
    T = _mod_eltype(var, mod)
    return FileArray{S,T,N}(filename, size(var); eachchunk, haschunks, mod, kw...)
end

# FileArray has S, T and N parameters not recoverable from fields
ConstructionBase.constructorof(::Type{<:FileArray{S,T,N}}) where {S,T,N} = FileArray{S,T,N}

filename(A::FileArray) = A.filename
mod(A::FileArray) = A.mod
DD.name(A::FileArray) = A.name
Base.size(A::FileArray) = A.size
DA.eachchunk(A::FileArray) = A.eachchunk
DA.haschunks(A::FileArray) = A.haschunks

# Run function `f` on the result of _open for the file type
function Base.open(f::Function, A::FileArray{S}; write=A.write, kw...) where S
    _open(f, S(), filename(A); name=name(A), group=A.group, write, mod=mod(A), kw...)
end

function DA.readblock!(A::FileArray, dst, r::AbstractUnitRange...)
    open(A) do O
        if isdisk(O)
            DA.readblock!(O, dst, r...)
        else
            dst[r...] .= view(parent(O), r...)
        end
    end
end
function DA.writeblock!(A::FileArray, src, r::AbstractUnitRange...)
    open(A; write=A.write) do O
        if isdisk(A)
            DA.writeblock!(O, src, r...)
        else
            parent(O)[r...] .= src
        end
    end
end


"""
    RasterDiskArray <: DiskArrays.AbstractDiskArray

A basic DiskArrays.jl wrapper for objects that don't have one defined yet.
When we `open` a `FileArray` it is replaced with a `RasterDiskArray`.
"""
struct RasterDiskArray{S,T,N,V,EC,HC,A} <: DiskArrays.AbstractDiskArray{T,N}
    var::V
    eachchunk::EC
    haschunks::HC
    attrib::A
end
function RasterDiskArray{S}(
    var::V, eachchunk::EC=DA.eachchunk(var), haschunks::HC=DA.haschunks(var), attrib::A=nothing
) where {S,V,EC,HC,A}
    T = eltype(var)
    N = ndims(var)
    RasterDiskArray{S,T,N,V,EC,HC,A}(var, eachchunk, haschunks, attrib)
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
MissingDiskArray{MT,N}(var::V) where {MT,N,V} = MissingDiskArray{MT,N,V}(var)

struct MissingDiskArrayConstructor{T,N} end
(::MissingDiskArrayConstructor{T,N})(var) where {T,N} = MissingDiskArray{T,N}(var)

ConstructionBase.constructorof(::Type{MissingDiskArray{T,N,V}}) where {T,N,V} =
    MissingDiskArrayConstructor{T,N}()

Base.parent(A::MissingDiskArray) = A.var
Base.size(A::MissingDiskArray) = size(parent(A))

DA.haschunks(A::MissingDiskArray) = DA.haschunks(parent(A))
DA.eachchunk(A::MissingDiskArray) = DA.eachchunk(parent(A))
DA.readblock!(A::MissingDiskArray, aout, r::AbstractUnitRange...) = DA.readblock!(parent(A), aout, r...)
DA.writeblock!(A::MissingDiskArray, v, r::AbstractUnitRange...) = DA.writeblock!(parent(A), v, r...)
