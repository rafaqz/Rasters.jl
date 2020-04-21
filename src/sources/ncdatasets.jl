using .NCDatasets

export NCDarray, NCDstack, NCDstackMetadata, NCDarrayMetadata, NCDdimMetadata

struct NCDstackMetadata{K,V} <: ArrayMetadata{K,V}
    val::Dict{K,V}
end

struct NCDarrayMetadata{K,V} <: ArrayMetadata{K,V}
    val::Dict{K,V}
end

struct NCDdimMetadata{K,V} <: DimMetadata{K,V}
    val::Dict{K,V}
end


# Utils ########################################################################

ncread(f, path::String) = NCDatasets.Dataset(f, path)

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

# Add a var array to a dataset before writing it.
ncwritevar!(dataset, A::AbstractGeoArray{T,N}) where {T,N} = begin
    A = reorderindex(A, Forward()) |>
        a -> reorderrelation(a, Forward())
    if ismissing(missingval(A))
        # TODO default _FillValue for Int?
        fillvalue = get(metadata(A), "_FillValue", NaN)
        A = replace_missing(A, convert(T, fillvalue))
    end
    # Define required dims
    for dim in dims(A)
        key = lowercase(name(dim))
        haskey(dataset.dim, key) && continue
        index = [val(dim)...]
        md = metadata(dim)
        attribvec = [] #md isa Nothing ? [] : [val(md)...]
        defDim(dataset, key, length(index))
        println("writing key: ", key, " of type: ", eltype(index))
        defVar(dataset, key, index, (key,); attrib=attribvec)
    end
    # TODO actually convert the metadata type
    attrib = if metadata isa NCDarrayMetadata
        deepcopy(val(metadata(A)))
    else
        Dict()
    end
    # Remove stack metdata if it is attached
    pop!(attrib, "dataset", nothing)
    # Set missing value
    if !ismissing(missingval(A))
        attrib["_FillValue"] = convert(T, missingval(A))
    end
    key = name(A)
    println("writing key: ", key, " of type: ", T)
    dimnames = lowercase.(name.(dims(A)))
    attribvec = [attrib...]
    var = defVar(dataset, key, eltype(A), dimnames; attrib=attribvec)

    if N == 0
        var[] = data(A)
    elseif N == 1
        var[:] = data(A)
    elseif N == 2
        var[:,:] = data(A)
    elseif N == 3
        var[:,:,:] = data(A)
    elseif N == 4
        var[:,:,:,:] = data(A)
    elseif N == 5
        var[:,:,:,:,:] = data(A)
    elseif N == 6
        var[:,:,:,:,:,:] = data(A)
    elseif N == 7
        var[:,:,:,:,:,:,:] = data(A)
    end

end

# CF standards don't enforce dimension names.
# But these are common, and should take care of most dims.
const dimmap = Dict("lat" => Lat,
                    "latitude" => Lat,
                    "lon" => Lon,
                    "long" => Lon,
                    "longitude" => Lon,
                    "time" => Ti,
                    "lev" => Vert,
                    "level" => Vert,
                    "vertical" => Vert,
                    "x" => X,
                    "y" => Y,
                    "z" => Z,
                   )

# Array ########################################################################
"""
    NCDarray(filename::AbstractString; name="", refdims=())

Create an array from a path to a netcdf file. The first non-dimension
layer of the file will be used as the array.

## Arguments
- `filename`: `String` pointing to a netcdf file.

## Keyword arguments
- `name`: Name for the array. Will use array key if not supplied.
- `refdims`: Add dimension position array was sliced from. Mostly used programatically.
  loading the array. Can save on disk load time for large files.
"""
struct NCDarray{T,N,A,D<:Tuple,R<:Tuple,Na<:AbstractString,Me,Mi,S,K
               } <: DiskGeoArray{T,N,D,LazyArray{T,N}}
    filename::A
    dims::D
    refdims::R
    name::Na
    metadata::Me
    missingval::Mi
    size::S
    key::K
end
NCDarray(filename::AbstractString, key...; kwargs...) =
    ncread(dataset -> NCDarray(dataset, filename, key...; kwargs...), filename)
NCDarray(dataset::NCDatasets.Dataset, filename, key=nothing;
         refdims=(),
         dims=nothing,
         name=nothing,
         metadata=nothing,
        ) = begin
    keys_ = nondimkeys(dataset)
    key = key isa Nothing || !(string(key) in keys_) ? first(keys_) : string(key)
    var = dataset[key]
    dims = dims isa Nothing ? GeoData.dims(dataset, key) : dims
    name = name isa Nothing ? string(key) : name
    metadata_ = metadata isa Nothing ? GeoData.metadata(var, GeoData.metadata(dataset)) : metadata
    missingval = missing
    size_ = map(length, dims)
    T = eltype(var)
    N = length(dims)

    NCDarray{T,N,typeof.((filename,dims,refdims,name,metadata_,missingval,size_,key))...
       }(filename, dims, refdims, name, metadata_, missingval, size_, key)
end

key(A::NCDarray) = A.key

# AbstractGeoArray methods

withsourcedata(f, A::NCDarray) =
    ncread(dataset -> f(dataset[string(key(A))]), filename(A))

# Base methods

Base.size(A::NCDarray) = A.size

"""
    Base.write(filename::AbstractString, ::Type{NCDarray}, s::AbstractGeoArray)

Write an NCDarray to a netcdf file using NCDatasets.jl
"""
Base.write(filename::AbstractString, ::Type{NCDarray}, A::AbstractGeoArray) = begin
    # Remove the dataset metadata
    stackmd = pop!(deepcopy(val(metadata(A))), "dataset", Dict())
    dataset = NCDatasets.Dataset(filename, "c"; attrib=stackmd)
    try
        ncwritevar!(dataset, A)
    finally
        close(dataset)
    end
end

# Stack ########################################################################

struct NCDstack{T,R,W,M,K} <: DiskGeoStack{T}
    filename::T
    refdims::R
    window::W
    metadata::M
    kwargs::K
end

"""
    NCDstack(filenames; refdims=(), window=(), metadata=nothing)

A lazy GeoStack that loads netcdf files using NCDatasets.jl

Create a stack from a list of filenames.

# Arguments
-`filenames`: `Vector` of `String` paths to netcdf files.

# Keyword arguments
- `refdims`: Add dimension position array was sliced from. Mostly used programatically.
- `window`: can be a tuple of Dimensions, selectors or regular indices.
- `metadata`: Add additional metadata as a `Dict`.
- `keys`: Keys for the layer in each file in filenames. If these do not match a layer
  the first layer will be used. This is also the default.

# Examples
```julia
multifile_stack = NCDstack([path1, path2, path3, path4])
```
"""
NCDstack(filenames::Union{Tuple,Vector}; refdims=(), window=(), metadata=nothing,
         keys=Tuple(Symbol.((ncread(ds -> first(nondimkeys(ds)), fp) for fp in filenames))),
         kwargs...) =
    GeoStack(NamedTuple{keys}(filenames), refdims, window, metadata, childtype=NCDarray, kwargs)

"""
    NCDstack(filename; refdims=(), window=(), metadata=nothing)

A lazy GeoStack that loads netcdf files using NCDatasets.jl

Create a stack from the filename of a netcdf file.

# Arguments
-`filename`: `String` path to a netcdf file.

# Keyword arguments
- `refdims`: Add dimension position array was sliced from. Mostly used programatically.
- `window`: can be a tuple of Dimensions, selectors or regular indices.
- `metadata`: Add additional metadata as a `Dict`.

# Examples
```julia
stack = NCDstack(filename; window=(Lat(Between(20, 40),))
stack[:soil_temperature]
```
"""
NCDstack(filename::AbstractString; refdims=(), window=(), 
         metadata=ncread(metadata, filename), kwargs...) =
    NCDstack(filename, refdims, window, metadata, kwargs)

childtype(::NCDstack) = NCDarray

# AbstractGeoStack methods

withsource(f, ::Type{NCDarray}, path::AbstractString, key=nothing) = ncread(f, path)
withsourcedata(f, ::Type{NCDarray}, path::AbstractString, key) =
    ncread(d -> f(d[string(key)]), path)

# Override the default to get the dims of the specific key
dims(::NCDstack, dataset, key::Key) = dims(dataset, key)

missingval(stack::NCDstack) = missing

# Base methods

Base.keys(stack::NCDstack{<:AbstractString}) =
    cleankeys(ncread(nondimkeys, getsource(stack)))

"""
    Base.write(filename::AbstractString, ::Type{NCDstack}, s::AbstractGeoStack)

Write an NCDstack to a single netcdf file, using NCDatasets.jl.

Currently `Dimension` metadata is not handled, and array metadata from other
array types is ignored.
"""
Base.write(filename::AbstractString, ::Type{NCDstack}, s::AbstractGeoStack) = begin
    dataset = NCDatasets.Dataset(filename, "c"; attrib=val(metadata(s)))
    try
        map(key -> ncwritevar!(dataset, s[key]), keys(s))
    finally
        close(dataset)
    end
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
            dimtype = get(dimmap, dimname, Dim{Symbol(dimname)})
            # Order: data is always forwards, we check the index order
            order = dvar[end] > dvar[1] ? Ordered(Forward(), Forward(), Forward()) :
                                          Ordered(Reverse(), Forward(), Reverse())

            # Assume the locus is at the center of the cell if boundaries aren't provided.
            # http://cfconventions.org/cf-conventions/cf-conventions.html#cell-boundaries

            if eltype(dvar) <: Number
                beginhalfcell = abs((dvar[2] - dvar[1]) * 0.5)
                endhalfcell = abs((dvar[end] - dvar[end-1]) * 0.5)
                bounds = if length(dvar) > 1
                    if isrev(indexorder(order))
                        dvar[end] - endhalfcell, dvar[1] + beginhalfcell
                    else
                        dvar[1] - beginhalfcell, dvar[end] + endhalfcell
                    end
                else
                    dvar[1], dvar[1]
                end
                locus = (dimtype <: TimeDim) ? Start() : Center()
                mode = Sampled(order, Irregular(bounds), Intervals(locus))
            elseif eltype(dvar) <: Dates.AbstractTime
                locus = Start()
                bounds = if length(dvar) > 1
                    if isrev(indexorder(order))
                        dvar[end], dvar[1] + (dvar[1] - dvar[2])
                    else
                        dvar[1], dvar[end] + (dvar[end] - dvar[end - 1])
                    end
                else
                    dvar[1], dvar[1]
                end
                mode = Sampled(order, Irregular(bounds), Intervals(locus))
            else
                mode = Sampled(order, Irregular(), Points())
            end

            meta = metadata(dvar)
            # Add the dim containing the dimension var array
            push!(dims, dimtype(dvar[:], mode, meta))
        else
            # The var doesn't exist. Maybe its `complex` or some other marker,
            # so make it a custom `Dim` with `NoIndex`
            push!(dims, Dim{Symbol(dimname)}(1:size(v, i), NoIndex(), nothing))
        end
    end
    (dims...,)
end

metadata(dataset::NCDatasets.Dataset) = NCDstackMetadata(Dict{String,Any}(dataset.attrib))
metadata(dataset::NCDatasets.Dataset, key::Key) = metadata(dataset[string(key)])
metadata(var::NCDatasets.CFVariable) = NCDarrayMetadata(Dict{String,Any}(var.attrib))
metadata(var::NCDatasets.CFVariable, stackmetadata::NCDstackMetadata) = begin
    md = metadata(var)
    md["dataset"] = stackmetadata
    md
end

missingval(var::NCDatasets.CFVariable{<:Union{Missing}}) = missing

# crs(dataset::NCDatasets.Dataset)
# crs(var::NCDatasets.CFVariable)
