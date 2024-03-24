"""
    FileStack{S,K}

    FileStack{S,K}(filename, types, sizes, eachchunk, haschunks, write)

A wrapper object that holds file pointer and size/chunking
metadata for a multi-layered stack stored in a single file, 
typically netcdf or hdf5.

`S` is a backend type like `NCDsource`, and `K` is a tuple of `Symbol` keys.
"""
struct FileStack{S,K,F<:AbstractString,T,SZ,EC,HC}
    filename::F
    types::T
    sizes::SZ
    eachchunk::EC
    haschunks::HC
    write::Bool
end
function FileStack{S,K}(
    filename::F, types::T, sizes::SZ, eachchunk::EC, haschunks::HC, write::Bool
   ) where {S,K,F,T,SZ,EC,HC}
    FileStack{S,K,F,T,SZ,EC,HC}(filename, types, sizes, eachchunk, haschunks, write)
end

# FileStack has `S` parameter that is not recoverable from fields.
ConstructionBase.constructorof(::Type{<:FileStack{S,K}}) where {S,K} = FileStack{S,K} 

filename(fs::FileStack) = fs.filename
Base.keys(fs::FileStack{<:Any,K}) where K = K
Base.values(fs::FileStack{<:Any}) = (fs[k] for k in keys(fs))
# Indexing FileStack returns a FileArray, 
# referencing a specific key in the same file.
function Base.getindex(fs::FileStack{S,K}, key::Symbol) where {S,K}
    is = NamedTuple{K}(ntuple(identity, length(K)))
    i = is[key]
    size = fs.sizes[i]
    eachchunk = fs.eachchunk[i]
    haschunks = fs.haschunks[i]
    T = fs.types[i]
    N = length(size)
    A = FileArray{S,T,N}(filename(fs), size, key, eachchunk, haschunks, fs.write)
end
