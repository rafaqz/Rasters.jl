"""
    FileStack{X,K}

    FileStack{X,K}(filename, types, sizes, eachchunk, haschunks, write)

A wrapper object that holds file pointer and size/chunking
metadata for a multi-layered stack stored in a single file, 
typically netcdf or hdf5.

`X` is a backend singleton like `GDALfile`, and `K` is a tuple
of `Symbol` keys.
"""
struct FileStack{X,K,F<:AbstractString,T<:NamedTuple{K},S<:NamedTuple{K},EC,HC}
    filename::F
    types::T
    sizes::S
    eachchunk::EC
    haschunks::HC
    write::Bool
end
function FileStack{X}(
    filename::F, types::T, sizes::S, eachchunk::EC, haschunks::HC, write::Bool
   ) where {X,K,F,T<:NamedTuple{K},S<:NamedTuple{K},EC,HC}
    FileStack{X,K,F,T,S,EC,HC}(filename, types, sizes, eachchunk, haschunks, write)
end

# FileStack has `X` parameter that is not recoverable from fields.
ConstructionBase.constructorof(::Type{<:FileStack{X}}) where {X} = FileStack{X} 

filename(fs::FileStack) = fs.filename
Base.keys(fs::FileStack{<:Any,K}) where K = K
Base.values(fs::FileStack{<:Any}) = (fs[k] for k in keys(fs))
# Indexing FileStack returns a FileArray, 
# referencing a specific key in the same file.
function Base.getindex(fs::FileStack{X}, key::Symbol) where X
    size = fs.sizes[key]
    eachchunk = fs.eachchunk[key]
    haschunks = fs.haschunks[key]
    T = fs.types[key]
    N = length(size)
    A = FileArray{X,T,N}(filename(fs), size, key, eachchunk, haschunks, fs.write)
end
