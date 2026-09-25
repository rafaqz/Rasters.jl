"""
    FileStack{S,Na}

    FileStack{S,Na}(filename, types, sizes, eachchunk, haschunks, write)

A wrapper object that holds a filepath string and size/chunking
metadata for a multi-layered stack stored in a single file. 
such as zarr or netcdf.

`S` is a backend type like `NCDsource`, and `Na` is a tuple of `Symbol` keys.
"""
struct FileStack{S,Na,T,SZ,G<:Union{AbstractString,Symbol,Nothing},EC,HC,M,OK<:NamedTuple}
    filename::String
    sizes::SZ
    group::G
    eachchunk::EC
    haschunks::HC
    mods::M
    write::Bool
    open_kw::OK
end
function FileStack{S,Na,T}(
    filename::AbstractString, sizes::SZ, group::G, eachchunk::EC, haschunks::HC, mods::M, write::Bool,
    open_kw::OK=NamedTuple(),
) where {S,Na,T,SZ,G,EC,M,HC,OK}
    FileStack{S,Na,T,SZ,G,EC,HC,M,OK}(String(filename), sizes, group, eachchunk, haschunks, mods, write, open_kw)
end
function FileStack{source}(ds::AbstractDataset, filename::AbstractString;
    write::Bool=false,
    group=nokw,
    name::NTuple{N,Symbol},
    vars,
    mods,
    open_kw=NamedTuple(),
) where {source,N}
    T = NamedTuple{name,Tuple{map(_mod_eltype, vars, mods)...}}
    layersizes = map(size, vars)
    eachchunk = map(DiskArrays.eachchunk, vars)
    haschunks = map(DiskArrays.haschunks, vars)
    group = isnokw(group) ? nothing : group
    return FileStack{source,name,T}(filename, layersizes, group, eachchunk, haschunks, mods, write, open_kw)
end

# FileStack has `S,Na,T` parameters that are not recoverable from fields.
ConstructionBase.constructorof(::Type{<:FileStack{S,Na,T}}) where {S,Na,T} = FileStack{S,Na,T} 

filename(fs::FileStack) = fs.filename
mods(fs::FileStack) = fs.mods

DD.name(::FileStack{<:Any,Na}) where Na = Na
DD.data_eltype(::FileStack{<:Any,<:Any,T}) where T = T

Base.eltype(::FileStack{<:Any,<:Any,T}) where T = T
Base.keys(::FileStack{<:Any,Na}) where Na = Na
Base.values(fs::FileStack{<:Any,Na}) where Na = (fs[n] for n in Na)

# Indexing FileStack returns a FileArray, 
# referencing a specific name in the same file.
function Base.getindex(fs::FileStack{S,Na,T}, name::Symbol) where {S,Na,T}
    is = NamedTuple{Na}(ntuple(identity, length(Na)))
    i = is[name]
    size = fs.sizes[i]
    eachchunk = fs.eachchunk[i]
    haschunks = fs.haschunks[i]
    mod = fs.mods[i]
    N = length(size)
    return FileArray{S,_itype(T, i),N}(filename(fs), size, name, fs.group, eachchunk, haschunks, mod, fs.write, fs.open_kw)
end

@inline _itype(::Type{<:NamedTuple{<:Any,T}}, i) where T = T.parameters[i]

