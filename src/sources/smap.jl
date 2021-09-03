using HDF5

export SMAPseries, SMAPstack, smapseries

const SMAPMISSING = -9999.0f0
const SMAPGEODATA = "Geophysical_Data"
const SMAPCRS = ProjString("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
const SMAPSIZE = (3856, 1624)
const SMAPDIMTYPES = (X, Y) 

# Dataset wrapper ###############################################################
# Becauase SMAP is just one of many HDF5 formats,
# we wrap it in SMAPhdf5 and SMAPvar wrappers

struct SMAPhdf5{T}
    ds::T
end

missingval(ds::SMAPhdf5) = SMAPMISSING
layermissingval(ds::SMAPhdf5) = SMAPMISSING
layerkeys(ds::SMAPhdf5) = keys(ds)
filekey(ds::SMAPhdf5, key::Nothing) = first(keys(ds))

function DD.dims(wrapper::SMAPhdf5)
    dataset = parent(wrapper)
    proj = read(HDF5.attributes(HDF5.root(dataset)["EASE2_global_projection"]), "grid_mapping_name")
    if proj == "lambert_cylindrical_equal_area"
        # There are matrices for lookup but all rows/colums are identical.
        # For performance and simplicity we just take a vector slice for each dim.
        extent = HDF5.attributes(HDF5.root(dataset)["Metadata/Extent"])
        lonbounds = read(extent["westBoundLongitude"]), read(extent["eastBoundLongitude"])
        latbounds = read(extent["southBoundLatitude"]), read(extent["northBoundLatitude"])
        latvec = HDF5.root(dataset)["cell_lat"][1, :]
        lonvec = HDF5.root(dataset)["cell_lon"][:, 1]
        lonmode = Mapped(
           order=Ordered(),
           span=Irregular(lonbounds),
           sampling=Intervals(Center()),
           crs=SMAPCRS,
           mappedcrs=EPSG(4326)
        )
        latmode = Mapped(
            order=Ordered(ReverseIndex(), ReverseArray(), ForwardRelation()),
            span=Irregular(latbounds),
            sampling=Intervals(Center()),
            crs=SMAPCRS,
            mappedcrs=EPSG(4326),
        )
        (X(lonvec; mode=lonmode), Y(latvec; mode=latmode))
    else
        error("projection $proj not supported")
    end
end

# TODO actually add metadata to the dict
DD.metadata(wrapper::SMAPhdf5) = Metadata{SMAPfile}(Dict())

function DD.layerdims(ds::SMAPhdf5)
    keys = cleankeys(layerkeys(ds))
    # All dims are the same
    NamedTuple{keys}(map(_ -> SMAPDIMTYPES, keys))
end

function DD.layermetadata(ds::SMAPhdf5)
    keys = cleankeys(layerkeys(ds))
    NamedTuple{keys}(map(_ -> DD.metadata(ds), keys))
end

Base.keys(ds::SMAPhdf5) = cleankeys(keys(parent(ds)[SMAPGEODATA]))
Base.parent(wrapper::SMAPhdf5) = wrapper.ds
Base.getindex(wrapper::SMAPhdf5, key) = SMAPvar(wrapper.ds[_smappath(key)])

_smappath(key::Key) = SMAPGEODATA * "/" * string(key)

struct SMAPvar{T}
    ds::T
end

Base.parent(wrapper::SMAPvar) = wrapper.ds
Base.eltype(wrapper::SMAPvar) = eltype(parent(wrapper))
Base.size(wrapper::SMAPvar) = SMAPSIZE
Base.ndims(wrapper::SMAPvar) = length(SMAPSIZE)
Base.getindex(wrapper::SMAPvar, args...) = getindex(parent(wrapper), args...)
Base.setindex!(wrapper::SMAPvar, args...) = setindex!(parent(wrapper), args...)

DA.eachchunk(var::SMAPvar) = DA.GridChunks(var, size(var))
DA.haschunks(var::SMAPvar) = DA.Unchunked()

# GeoArray ######################################################################

function FileArray(ds::SMAPhdf5, filename::AbstractString; key, kw...)
    FileArray(ds[key], filename; key, kw...)
end
function FileArray(var::SMAPvar, filename::AbstractString; key, kw...)
    T = eltype(var)
    N = ndims(var)
    eachchunk = DA.eachchunk(var)
    haschunks = DA.haschunks(var)
    FileArray{SMAPfile,T,N}(filename, SMAPSIZE; key, eachchunk, haschunks, kw...)
end

function Base.open(f::Function, A::FileArray{SMAPfile}; kw...)
    _open(SMAPfile, filename(A); key=key(A), kw...) do var
        f(GeoDiskArray{SMAPfile}(var)) 
    end
end
    
DA.writeblock!(A::GeoDiskArray{SMAPfile}, v, r::AbstractUnitRange...) = A[r...] = v

haslayers(::Type{SMAPfile}) = true

# Stack ########################################################################

@deprecate SMAPstack(args...; kw...) GeoStack(args...; source=SMAPfile, kw...)

function FileStack{SMAPfile}(ds::SMAPhdf5, filename::AbstractString; write=false, keys)
    keys = map(Symbol, keys isa Nothing ? layerkeys(ds) : keys) |> Tuple
    type_size_ec_hc = map(keys) do key
        var = GeoDiskArray{SMAPfile}(ds[key])
        eltype(var), size(var), DA.eachchunk(var), DA.haschunks(var)
    end
    layertypes = NamedTuple{keys}(map(x->x[1], type_size_ec_hc))
    layersizes = NamedTuple{keys}(map(x->x[2], type_size_ec_hc))
    eachchunk = NamedTuple{keys}(map(x->x[3], type_size_ec_hc))
    haschunks = NamedTuple{keys}(map(x->x[4], type_size_ec_hc))
    FileStack{SMAPfile}(filename, layertypes, layersizes, eachchunk, haschunks, write)
end

# Series #######################################################################

"""
    SMAPseries(filenames::AbstractString; kw...)
    SMAPseries(filenames::Vector{<:AbstractString}, dims=nothing; kw...)

[`GeoSeries`](@ref) loader for SMAP files and whole folders of files,
organised along the time dimension. Returns a [`GeoSeries`](@ref).

# Arguments

- `filenames`: A `String` path to a directory of SMAP files,
    or a vector of `String` paths to specific files.
- `dims`: `Tuple` containing `Ti` dimension for the series.
    Automatically generated form `filenames` unless passed in.

# Keywords

- `kw`: Passed to `GeoSeries`.
"""
function smapseries(dir::AbstractString; kw...)
    smapseries(joinpath.(dir, filter_ext(dir, ".h5")); kw...)
end
function smapseries(filenames::Vector{<:AbstractString}, dims=nothing; kw...)
    if dims isa Nothing
        usedpaths = String[]
        timeseries = []
        errors = []
        for filename in filenames
            try
                t = _smap_timefromfilename(filename)
                push!(timeseries, t)
                push!(usedpaths, filename)
            catch e
                push!(errors, e)
            end
        end
        # Use the first files time dim as a template, but join vals into an array of times.
        timedim = _smap_timedim(timeseries)
    else
        usedpaths = filenames
    end
    # Show errors after all files load, or you can't see them:
    if length(errors) > 0
        println("Some errors thrown during file load: ")
        println.(errors)
    end
    # Get the dims once for the whole series
    dims, metadata =_open(SMAPfile, first(filenames)) do ds
        DD.dims(ds), DD.metadata(ds)
    end
    GeoSeries(usedpaths, (timedim,); child=stack, dims, metadata, kw...)
end

@deprecate SMAPseries(args...; kw...) smapseries(args...; kw...)


# Utils ########################################################################

function _open(f, ::Type{SMAPfile}, filepath::AbstractString; key=nothing, kw...)
    if key isa Nothing
        h5open(filepath; kw...) do ds
            cleanreturn(f(SMAPhdf5(ds)))
        end
    else
        h5open(filepath) do ds
            cleanreturn(f(SMAPhdf5(ds)[key]))
        end
    end
end

function _smap_timefromfilename(filename::String)
    dateformat = DateFormat("yyyymmddTHHMMSS")
    dateregex = r"SMAP_L4_SM_gph_(\d+T\d+)_"
    datematch = match(dateregex, filename)
    if datematch !== nothing
        return DateTime(datematch.captures[1], dateformat)
    else
        error("Date/time not correctly formatted in path: $filename")
    end
end

_smap_timedim(t::DateTime) = _smap_timedim(t:Hour(3):t)
_smap_timedim(times::AbstractVector) =
    Ti(times, mode=Sampled(Ordered(), Regular(Hour(3)), Intervals(Start())))
