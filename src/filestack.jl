"""
    FileStack{S,Na}

    FileStack{S,Na}(filename, types, sizes, eachchunk, haschunks, write)

A wrapper object that holds file pointer and size/chunking
metadata for a multi-layered stack stored in a single file, 
typically netcdf or hdf5.

`S` is a backend type like `NCDsource`, and `Na` is a tuple of `Symbol` keys.
"""
struct FileStack{S,Na,T,SZ,G<:Union{AbstractString,Symbol,Nothing},EC,HC}
    filename::String
    types::T
    sizes::SZ
    group::G
    eachchunk::EC
    haschunks::HC
    write::Bool
end
function FileStack{S,Na}(
    filename, types::T, sizes::SZ, group::G, eachchunk::EC, haschunks::HC, write::Bool
   ) where {S,Na,T,SZ,G,EC,HC}
    FileStack{S,Na,T,SZ,G,EC,HC}(String(filename), types, sizes, group, eachchunk, haschunks, write)
end

# FileStack has `S` parameter that is not recoverable from fields.
ConstructionBase.constructorof(::Type{<:FileStack{S,Na}}) where {S,Na} = FileStack{S,Na} 

filename(fs::FileStack) = fs.filename
DD.name(fs::FileStack{<:Any,Na}) where Na = Na
Base.keys(fs::FileStack{<:Any,Na}) where Na = Na
Base.values(fs::FileStack{<:Any}) = (fs[k] for k in keys(fs))
# Indexing FileStack returns a FileArray, 
# referencing a specific name in the same file.
function Base.getindex(fs::FileStack{S,Na}, name::Symbol) where {S,Na}
    is = NamedTuple{Na}(ntuple(identity, length(Na)))
    i = is[name]
    size = fs.sizes[i]
    eachchunk = fs.eachchunk[i]
    haschunks = fs.haschunks[i]
    T = fs.types[i]
    N = length(size)
    return FileArray{S,T,N}(filename(fs), size, name, fs.group, eachchunk, haschunks, fs.write)
end
