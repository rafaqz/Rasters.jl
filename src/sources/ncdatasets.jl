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
    for (i, dimname) in enumerate(NCDatasets.dimnames(v))
        if haskey(ds, dimname)
            dvar = ds[dimname]
            # Find the matching dimension constructor. If its an unknown name use 
            # the generic Dim with the dim name as type parameter
            dimconstructor = get(dimmap, dimname, Dim{Symbol(dimname)})
            # Get the attrib metadata
            meta = Dict(metadata(dvar))
            # Add the dim containing the dimension var array 
            push!(dims, dimconstructor(dvar[:], meta))
        else
            # The var doesn't exist. Maybe its `complex` or some other marker
            # so just make it a Dim with that name and range matching the indices
            push!(dims, Dim{Symbol(dimname)}(1:size(v, i)))
        end
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
NCarray(ds::NCDatasets.Dataset, key=first(nondimkeys(ds)); refdims=()) = begin
    var = ds[string(key)]
    NCarray(Array(var), dims(ds, key), refdims, metadata(var), missingval(var))
end


struct NCstack{T} <: AbstractGeoStack{T}
    data::T
end
NCstack(filepaths::Union{Tuple,Vector}) = begin
    keys = Tuple(Symbol.((NCDatasets.Dataset(ds->first(nondimkeys(ds)), fp) for fp in filepaths)))
    NCstack(NamedTuple{keys}(filepaths))
end

run(f, stack::NCstack, path) = NCDatasets.Dataset(f, path)
data(stack::NCstack, ds, key::Key) = NCarray(ds, key)
data(stack::NCstack, ds, key::Key, I...) = NCarray(ds, key)[I...] 
dims(stack::NCstack, ds, key::Key) = dims(ds, key)
metadata(stack::NCstack) = nothing
missingval(stack::NCstack) = missing
refdims(stack::NCstack) = ()

@inline Base.keys(stack::AbstractGeoStack{<:AbstractString}) = run(nondimkeys, stack, source(stack))
@inline Base.keys(stack::AbstractGeoStack{<:NamedTuple}) = keys(parent(stack))
Base.copy!(dst::AbstractArray, src::NCstack, key) = 
    run(ds -> copy!(dst, ds[string(key)]), src, source(src))

# utils

nondimkeys(ds) = begin
    dimkeys = keys(ds.dim)
    if "bnds" in dimkeys
        dimkeys = setdiff(dimkeys, ("bnds",))
        boundskeys = (k -> ds[k].attrib["bounds"]).(dimkeys)
        dimkeys = union(dimkeys, boundskeys)
    end
    setdiff(keys(ds), dimkeys)
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
