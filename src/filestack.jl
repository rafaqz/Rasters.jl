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
    sizes::SZ
    group::G
    eachchunk::EC
    haschunks::HC
    write::Bool
end
function FileStack{S,Na,T}(
    filename::AbstractString, sizes::SZ, group::G, eachchunk::EC, haschunks::HC, write::Bool
) where {S,Na,T,SZ,EC,HC}
    FileStack{S,Na,T,SZ,EC,HC}(String(filename), sizes, group, eachchunk, haschunks, write)
end

# FileStack has `S,Na,T` parameters that are not recoverable from fields.
ConstructionBase.constructorof(::Type{<:FileStack{S,K,T}}) where {S,Na,T} = FileStack{S,Na,T} 

filename(fs::FileStack) = fs.filename

DD.name(::FileStack{<:Any,Na}) where Na = Na
DD.data_eltype(::FileStack{<:Any,<:Any,T}) where T = T

Base.eltype(::FileStack{<:Any,<:Any,T}) where T = T
Base.keys(::FileStack{<:Any,Na}) where Na = Na
Base.values(fs::FileStack{<:Any,Na}) where Na = (fs[n] for n in Na)

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
