using NCDatasets

export NCDarray, NCDstack, NCDmetadata, NCDdimMetadata

struct NCDmetadata{M} <: AbstractArrayMetadata
    val::M
end

struct NCDdimMetadata{M} <: AbstractDimMetadata
    val::M
end

# Array ########################################################################

struct NCDarray{T,N,A,D<:Tuple,R<:Tuple,Me,Mi,Na,W,S} <: AbstractGeoArray{T,N,D}
    filename::A
    dims::D
    refdims::R
    metadata::Me
    missingval::Mi
    name::Na
    window::W
    size::S
end

NCDarray(filename::AbstractString; kwargs...) =
    ncapply(dataset -> NCDarray(dataset, filename; kwargs...), filename)
NCDarray(dataset::NCDatasets.Dataset, filename;
         dims=dims(dataset),
         refdims=(),
         name=first(nondimkeys(dataset)),
         metadata = metadata(dataset),
         window=()) = begin
    var = dataset[name]
    if window == ()
        sze = size(var)
    else
        sze = windowsize(window)
        dims, refdims = slicedims(dims, refdims, window)
    end
    missingval = missing
    T = eltype(var)
    N = ndims(var)
    NCDarray{T,N,typeof.((filename,dims,refdims,metadata,missingval,name,window,sze))...
       }(filename, dims, refdims, metadata, missingval, name, window, sze)
end

Base.size(A::NCDarray) = A.size
Base.parent(A::NCDarray) = ncapply(A.filename) do dataset
    ncread(dataset[name(A)], windoworempty(A)...)
end
Base.getindex(A::NCDarray, I::Vararg{<:Union{<:Integer,<:AbstractArray}}) = begin
    I = applywindow(A, I)
    rebuildsliced(A, ncapply(dataset -> ncread(dataset[name(A)], I...), A.filename), I)
end
Base.getindex(A::NCDarray, I::Vararg{<:Integer}) = begin
    I = applywindow(A, I)
    ncapply(dataset -> ncread(dataset[name(A)], I...), A.filename)
end


# Stack ########################################################################

struct NCDstack{T,D,R,W,M} <: DiskGeoStack{T}
    data::T
    dims::D
    refdims::R
    window::W
    metadata::M
end

"""
    NCDstack(filepaths::Union{Tuple,Vector}; dims=(), refdims=(), window=(), metadata=Nothing)

Create a stack from an array or tuple of paths to netcdf files. The first non-dimension
layer of each file will be used in the stack.

This constructor is intended for handling simple single-layer netcdfs.
"""
NCDstack(filepaths::Union{Tuple,Vector}; dims=ncapply(dims, first(filepaths)),
        refdims=(), window=(), metadata=Nothing) = begin
    keys = Tuple(Symbol.((ncapply(dataset->first(nondimkeys(dataset)), fp) for fp in filepaths)))
    NCDstack(NamedTuple{keys}(filepaths), dims, refdims, window, metadata)
end
NCDstack(data::String; dims=ncapply(dims, data),
        refdims=(), window=(), metadata=ncapply(metadata, data)) =
    NCDstack(data, dims, refdims, window, metadata)


safeapply(f, ::NCDstack, path) = ncapply(f, path)

data(s::NCDstack, dataset, key::Key, I...) =
    GeoArray(dataset[string(key)][I...], slicedims(dims(s, key), refdims(s), I)...,
             metadata(s), missingval(s), Symbol(key))
data(::NCDstack, dataset, key::Key, I::Vararg{Integer}) = dataset[string(key)][I...]
data(s::NCDstack, dataset, key::Key) =
    GeoArray(Array(dataset[string(key)]), dims(dataset, key), refdims(s),
             metadata(s), missingval(s), Symbol(key))

dims(::NCDstack, dataset, key::Key) = dims(dataset, key)
dims(::NCDstack, dataset, key::Key) = dims(dataset, key)

missingval(stack::NCDstack) = missing

Base.keys(stack::NCDstack{<:AbstractString}) =
    Tuple(Symbol.(safeapply(nondimkeys, stack, source(stack))))
Base.copy!(dst::AbstractArray, src::NCDstack, key) =
    safeapply(dataset -> copy!(dst, dataset[string(key)]), src, source(src))
Base.copy!(dst::AbstractGeoArray, src::NCDstack, key::Key) =
    safeapply(dataset -> copy!(parent(dst), dataset[string(key)]), src, source(src))

Base.write(filename::AbstractString, ::Type{NCDstack}, s::AbstractGeoStack) = begin
    dataset = NCDatasets.Dataset(filename, "c"; attrib=val(metadata(s)))
    try
        map(key -> addvar!(dataset, s[key]), keys(s))
    finally
        close(dataset)
    end
end

Base.write(filename::AbstractString, ::Type{NCDarray}, A::AbstractGeoArray) = begin
    dataset = NCDatasets.Dataset(filename, "c")
    try
        addvar!(dataset, A)
    finally
        close(dataset)
    end
end

addvar!(dataset, A) = begin
    A = forwardorder(A)
    if ismissing(missingval(A))
        fillvalue = get(metadata(A), "_FillValue", NaN)
        A = replace_missing(A, convert(eltype(A), fillvalue))
    end

    # Define required dims
    for dim in dims(A)
        key = lowercase(name(dim))
        haskey(dataset.dim, key) && continue
        index = [val(dim)...]
        defDim(dataset, key, length(index))
        # if isnothing(metadata(dim))
        println("writing key: ", string(key))
        defVar(dataset, key, index, (key,))
        # TODO this hangs.
        # if typeof(dim) <: DimensionalData.Time
            # defVar(dataset, key, index, (key,); attrib=[metadata(dim)...])
        # end
    end
    attrib = Dict()
    if !ismissing(missingval(A))
        attrib["_FillValue"] = convert(eltype(A), missingval(A))
    end
    # println("FillValue: ", attrib["_FillValue"])
    _name = string(name(A))
    if _name == "" _name = "Unnamed" end
    println("writing key: ", _name, " of type: ", eltype(parent(A)))
    var = defVar(dataset, _name, eltype(parent(A)), lowercase.(name.(dims(A)));
                 attrib=[attrib...])

    var[:] = parent(A)
end

# DimensionalData methods for NCDatasets types ###############################

dims(dataset::NCDatasets.Dataset) = dims(dataset, first(nondimkeys(dataset)))
dims(dataset::NCDatasets.Dataset, key::Key) = begin
    v = dataset[string(key)]
    dims = []
    for (i, dimname) in enumerate(NCDatasets.dimnames(v))
        if haskey(dataset, dimname)
            dvar = dataset[dimname]
            # Find the matching dimension constructor. If its an unknown name use
            # the generic Dim with the dim name as type parameter
            dimconstructor = get(dimmap, dimname, Dim{Symbol(dimname)})
            # Get the attrib metadata
            order = dvar[end] > dvar[1] ? Ordered(Forward(), Forward(), Forward()) : 
                                          Ordered(Reverse(), Reverse(), Forward())
            grid = AllignedGrid(order=order)
            meta = metadata(dvar)
            # Add the dim containing the dimension var array
            push!(dims, dimconstructor(dvar[:], grid, meta))
        else
            # The var doesn't exist. Maybe its `complex` or some other marker
            # so just make it a Dim with that name and range matching the indices
            push!(dims, Dim{Symbol(dimname)}(1:size(v, i)))
        end
    end
    dims = formatdims(v, (dims...,))
end
metadata(dataset::NCDatasets.Dataset) = NCDmetadata(Dict(dataset.attrib))
metadata(dataset::NCDatasets.Dataset, key::Key) = metadata(dataset[string(key)])

metadata(var::NCDatasets.CFVariable) = NCDmetadata(Dict(var.attrib))
missingval(var::NCDatasets.CFVariable{<:Union{Missing}}) = missing


# Utils ########################################################################

# CF standards don't enforce dimension names.
# But these are common, and should take care of most dims.
const dimmap = Dict("lat" => Lat,
                    "latitude" => Lat,
                    "lon" => Lon,
                    "long" => Lon,
                    "longitude" => Lon,
                    "time" => Time,
                    "lev" => Vert,
                    "level" => Vert,
                    "vertical" => Vert,
                    "x" => X,
                    "y" => Y,
                    "z" => Z,
                   )

ncapply(f, path::String) = NCDatasets.Dataset(f, path)

ncread(var, I...) = var[I...]
ncread(var) = Array(var)


nondimkeys(dataset) = begin
    dimkeys = keys(dataset.dim)
    removekeys = if "bnds" in dimkeys
        dimkeys = setdiff(dimkeys, ("bnds",))
        boundskeys = map(k -> dataset[k].attrib["bounds"], dimkeys)
        union(dimkeys, boundskeys)
    else
        dimkeys
    end
    setdiff(keys(dataset), removekeys)
end
