# module NCDgeoData

using NCDatasets

export NCDarray, NCDstack

# CF standards don't enforce dimension names. 
# But these are common, and should take care most dims.
const dimmap = Dict("lat" => Lat, 
                    "latitude" => Lat, 
                    "lon" => Lon, 
                    "long" => Lon, 
                    "longitude" => Lon, 
                    "time" => Time, 
                    "lev" => Vert, 
                    "level" => Vert, 
                    "vertical" => Vert) 

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
            order = dvar[end] > dvar[1] ? Ordered(Forward(), Forward()) : Ordered(Reverse(), Reverse())
            grid = AllignedGrid(order=order)
            meta = Dict(metadata(dvar))
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
metadata(dataset::NCDatasets.Dataset) = Dict(dataset.attrib)
metadata(dataset::NCDatasets.Dataset, key::Key) = metadata(dataset[string(key)])

metadata(var::NCDatasets.CFVariable) = Dict(var.attrib)
missingval(var::NCDatasets.CFVariable{<:Union{Missing}}) = missing


# Array ########################################################################

@GeoArrayMixin struct NCDarray{A<:AbstractArray{T,N},W} <: AbstractGeoArray{T,N,D} 
    window::W
end

# TODO make this lazy?
NCDarray(path::AbstractString; refdims=(), window=()) = 
    ncapply(dataset -> NCDarray(dataset; refdims=refdims, window=window), path) 
NCDarray(dataset::NCDatasets.Dataset, key=first(nondimkeys(dataset)); 
        refdims=(), name=Symbol(key), window=()) = begin
    var = dataset[string(key)]
    NCDarray(Array(var), dims(dataset, key), refdims, metadata(var), missingval(var), name, window)
end


# Stack ########################################################################


@GeoStackMixin struct NCDstack{} <: AbstractGeoStack{T} end

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


# Utils ########################################################################

ncapply(f, path) = NCDatasets.Dataset(f, path)

nondimkeys(dataset) = begin
    dimkeys = keys(dataset.dim)
    if "bnds" in dimkeys
        dimkeys = setdiff(dimkeys, ("bnds",))
        boundskeys = (k -> dataset[k].attrib["bounds"]).(dimkeys)
        dimkeys = union(dimkeys, boundskeys)
    end
    setdiff(keys(dataset), dimkeys)
end

# save(s::NCDstack, path) = begin
#     dataset = Dataset(path, "c")

#     for (key, val) in metadata(s)
#         dataset[key] = val 
#     end

#     for (key, layer) in s
#         dimstrings = (shortname(dim) for dim in dims(value))

#         defDim.(Ref(dataset), dimstrings, size.(val.(dims.(layer))))

#         # Define a variable
#         v = defVar(dataset, key, eltype(v), size(layer))
#         # TODO: add dims to variable

#         for (key, val) in metadata(layer)
#             metadata(v)[key] = val
#         end
#         v .= replace_missing(parent(layer), NaN)
#     end
#     close(dataset)
# end
