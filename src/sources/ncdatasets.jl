using .NCDatasets

export NCDarray, NCDstack, NCDmetadata, NCDdimMetadata

struct NCDmetadata{K,V} <: ArrayMetadata{K,V}
    val::Dict{K,V}
end

struct NCDdimMetadata{K,V} <: DimMetadata{K,V}
    val::Dict{K,V}
end

# Array ########################################################################

struct NCDarray{T,N,A,D<:Tuple,R<:Tuple,Me,Mi,Na,W,S} <: DiskGeoArray{T,N,D,LazyArray{T,N}}
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
        window = dims2indices(dims, window)
        sze = windowsize(window)
    end
    missingval = missing
    T = eltype(var)
    N = length(sze)
    NCDarray{T,N,typeof.((filename,dims,refdims,metadata,missingval,name,window,sze))...
       }(filename, dims, refdims, metadata, missingval, name, window, sze)
end


data(A::NCDarray) =
    ncapply(filename(A)) do dataset
        var = dataset[name(A)]
        _window = maybewindow2indices(var, dims(A), window(A))
        ncread(var, _window)
    end
filename(A::NCDarray) = A.filename
crs(A::NCDarray) = ncapply(crs, filename(A))

Base.size(A::NCDarray) = A.size
Base.getindex(A::NCDarray, I::Vararg{<:Union{<:Integer,<:AbstractArray}}) =
    ncapply(filename(A)) do dataset
        var = dataset[name(A)]
        _window = maybewindow2indices(var, dims(A), window(A))
        # Slice for both window and indices
        _dims, _refdims = slicedims(slicedims(dims(A), refdims(A), _window)..., I)
        data = ncread(var, _window, I...)
        rebuild(A, data, _dims, _refdims)
    end
Base.getindex(A::NCDarray, I::Vararg{<:Integer}) =
    ncapply(filename(A)) do dataset
        var = dataset[name(A)]
        _window = maybewindow2indices(var, dims(A), window(A))
        ncread(var, _window, I...)
    end

Base.write(filename::AbstractString, ::Type{NCDarray}, A::AbstractGeoArray) = begin
    dataset = NCDatasets.Dataset(filename, "c")
    try
        ncaddvar!(dataset, A)
    finally
        close(dataset)
    end
end

# Stack ########################################################################

struct NCDstack{T,R,W,M} <: DiskGeoStack{T}
    filename::T
    refdims::R
    window::W
    metadata::M
end

"""
    NCDstack(filepaths::Union{Tuple,Vector}; refdims=(), window=(), metadata=Nothing)

Create a stack from an array or tuple of paths to netcdf files. The first non-dimension
layer of each file will be used in the stack.

This constructor is intended for handling simple single-layer netcdfs.
"""
NCDstack(filepaths::Union{Tuple,Vector}; refdims=(), window=(), metadata=Nothing,
         keys=Tuple(Symbol.((ncapply(ds -> first(nondimkeys(ds)), fp) for fp in filepaths)))) =
    NCDstack(NamedTuple{keys}(filepaths), refdims, window, metadata)
NCDstack(filename::String; refdims=(), window=(), metadata=ncapply(metadata, filename)) =
    NCDstack(filename, refdims, window, metadata)

safeapply(f, ::NCDstack, path) = ncapply(f, path)

@inline Base.getindex(s::NCDstack, key::Key, i1::Integer, I::Integer...) =
    ncapply(filename(s, key)) do dataset
        key = string(key)
        var = dataset[key]
        _window = maybewindow2indices(var, dims(dataset, key), window(s))
        ncread(var, _window, i1, I...)
    end
@inline Base.getindex(s::NCDstack, key::Key, I::Union{Colon,Integer,AbstractArray}...) =
    ncapply(filename(s, key)) do dataset
        key = string(key)
        var = dataset[key]
        _dims = dims(dataset, key)
        _window = maybewindow2indices(var, _dims, window(s))
        _dims, _refdims = slicedims(slicedims(_dims, refdims(s), _window)..., I)
        A = ncread(var, _window, I...)
        GeoArray(A, _dims, _refdims, metadata(s), missingval(s), key)
    end

dims(::NCDstack, dataset, key::Key) = dims(dataset, key)
dims(::NCDstack, dataset, key::Key) = dims(dataset, key)

missingval(stack::NCDstack) = missing

Base.keys(stack::NCDstack{<:AbstractString}) =
    Tuple(Symbol.(safeapply(nondimkeys, stack, source(stack))))

Base.copy!(dst::AbstractGeoArray, src::NCDstack, key::Key) =
    copy!(data(dst), src, key)
Base.copy!(dst::AbstractArray, src::NCDstack, key) =
    ncapply(filename(src)) do dataset
        key = string(key)
        var = dataset[key]
        _window = maybewindow2indices(var, dims(dataset, key), window(src))
        copy!(dst, readwindowed(var, _window))
    end

Base.write(filename::AbstractString, ::Type{NCDstack}, s::AbstractGeoStack) = begin
    dataset = NCDatasets.Dataset(filename, "c"; attrib=val(metadata(s)))
    try
        map(key -> ncaddvar!(dataset, s[key]), keys(s))
    finally
        close(dataset)
    end
end

ncaddvar!(dataset, A) = begin
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
        println("writing key: ", string(key))
        defVar(dataset, key, index, (key,))
    end
    attrib = Dict()
    if !ismissing(missingval(A))
        attrib["_FillValue"] = convert(eltype(A), missingval(A))
    end
    # println("FillValue: ", attrib["_FillValue"])
    _name = string(name(A))
    if _name == "" _name = "Unnamed" end
    println("writing key: ", _name, " of type: ", eltype(data(A)))
    var = defVar(dataset, _name, eltype(data(A)), lowercase.(name.(dims(A)));
                 attrib=[attrib...])

    var[:] = data(A)
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

            # Assume the locus is at the center of the cell if boundaries aren't provided.
            # http://cfconventions.org/cf-conventions/cf-conventions.html#cell-boundaries

            if eltype(dvar) isa Number
                beginhalfcell = abs((dvar[2] - dvar[1]) * 0.5)
                endhalfcell = abs((dvar[end] - dvar[end-1]) * 0.5)
                bounds = if isrev(indexorder(order))
                    dvar[end] - endhalfcell, dvar[1] + beginhalfcell
                else
                    dvar[1] - beginhalfcell, dvar[end] + endhalfcell
                end
                grid = BoundedGrid(order=order, locus=Center(), bounds=bounds)
            else
                grid = AlignedGrid(order=order)
            end


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

# crs(dataset::NCDatasets.Dataset)
# crs(var::NCDatasets.CFVariable) = NCDmetadata(Dict(var.attrib))


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

ncread(A, window::Tuple{}) = Array(A)
ncread(A, window::Tuple{}, I...) = A[I...]
ncread(A, window, I...) = A[Base.reindex(window, I)...]
ncread(A, window) = A[window...]


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
