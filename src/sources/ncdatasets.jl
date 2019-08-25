using NCDatasets

export NCarray, NCstack

# CF standards don't enforce dimension names. 
# These are common and should take care most common dims.
const dimmap = Dict("lat" => Lat, 
                    "latitude" => Lat, 
                    "lon" => Lon, 
                    "long" => Lon, 
                    "longitude" => Lon, 
                    "time" => Time, 
                    "lev" => Vert, 
                    "level" => Vert, 
                    "vertical" => Vert) 

# Methods for NCDatasets types

dims(ds::NCDatasets.Dataset) = dims(first(keys(ds)))
dims(ds::NCDatasets.Dataset, key::Key) = begin
    v = ds[string(key)]
    dims = []
    for dname in NCDatasets.dimnames(v)
        dvar = ds[dname]
        # Find the matching dimension constructor. If its an unknown name use 
        # the generic Dim with the dim name as type parameter
        dimconstructor = get(dimmap, dname, Dim{Symbol(dname)})
        # Get the attrib metadata
        meta = Dict(metadata(dvar))
        # Add the dim containing the dimension var array 
        push!(dims, dimconstructor(dvar[:], meta))
    end
    dims = formatdims(v, (dims...,))
end
metadata(ds::NCDatasets.Dataset) = ds.attrib
metadata(ds::NCDatasets.Dataset, key::Key) = metadta(ds[key])

metadata(var::NCDatasets.CFVariable) = var.attrib
missingval(var::NCDatasets.CFVariable{<:Union{Missing}}) = missing


@GeoArrayMixin struct NCarray{} <: AbstractGeoArray{T,N,D} end

NCarray(path::AbstractString; refdims=()) = 
    NCDatasets.Dataset(ds -> NCarray(ds; refdims=refdims), path) 
NCarray(ds::NCDatasets.Dataset, key=nondimkey(ds); refdims=()) = begin
    var = ds[string(key)]
    NCarray(Array(var), dims(ds, key), refdims, metadata(var), missingval(var))
end


struct NCstack{T} <: AbstractGeoStack{T}
    data::T
end
NCstack(filepaths::Union{Tuple,Vector}) = begin
    keys = Tuple(Symbol.((NCDatasets.Dataset(nondimkey, fp) for fp in filepaths)))
    NCstack(NamedTuple{keys}(filepaths))
end

run(f, stack::NCstack, path) = NCDatasets.Dataset(f, path)
data(stack::NCstack, ds, key::Key) = NCarray(ds, key)
data(stack::NCstack, ds, key::Key, I...) = NCarray(ds, key)[I...] 
dims(stack::NCstack, ds, key::Key) = dims(ds, key)
missingval(stack::NCstack, args...) = missing
refdims(stack::NCstack) = ()


# utils

nondimkey(ds) = begin
    dimkeys = keys(ds.dim)
    if "bnds" in dimkeys
        dimkeys = setdiff(dimkeys, ("bnds",))
        boundskeys = (k -> ds[k].attrib["bounds"]).(dimkeys)
        dimkeys = union(dimkeys, boundskeys)
    end
    nondimkeys = setdiff(keys(ds), dimkeys)
    key = nondimkeys[1]
end

# save(s::NCstack, path) = begin
#     ds = Dataset(path, "c")

#     for (key, val) in metadata(s)
#         ds[key] = val 
#     end

#     for (key, layer) in s
#         dimstrings = (shortname(dim) for dim in dims(value))

#         defDim.(Ref(ds), dimstrings, size.(val.(dims.(layer))))

#         # Define a variable
#         v = defVar(ds, key, eltype(v), size(layer))
#         # TODO: add dims to variable

#         for (key, val) in metadata(layer)
#             metadata(v)[key] = val
#         end
#         v .= replace_missing(parent(layer), NaN)
#     end
#     close(ds)
# end
