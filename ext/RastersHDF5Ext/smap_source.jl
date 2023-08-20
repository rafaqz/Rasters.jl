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

RA.missingval(ds::SMAPhdf5) = SMAPMISSING
RA.layerkeys(ds::SMAPhdf5) = keys(ds)
RA.filekey(ds::SMAPhdf5, key::Nothing) = first(keys(ds))

function DD.dims(wrapper::SMAPhdf5)
    dataset = parent(wrapper)
    proj = read(HDF5.attributes(HDF5.root(dataset)["EASE2_global_projection"]), "grid_mapping_name")
    if proj == "lambert_cylindrical_equal_area"
        # There are matrices for lookup but all rows/colums are identical.
        # For performance and simplicity we just take a vector slice for each dim.
        extent = HDF5.attributes(HDF5.root(dataset)["Metadata/Extent"])
        lonbounds = read(extent["westBoundLongitude"]), read(extent["eastBoundLongitude"])
        latbounds = read(extent["southBoundLatitude"]), read(extent["northBoundLatitude"])
        lonvec = HDF5.root(dataset)["cell_lon"][:, 1]
        latvec = HDF5.root(dataset)["cell_lat"][1, :]
        lonlookup = Mapped(lonvec;
           order=ForwardOrdered(),
           span=Irregular(lonbounds),
           sampling=Intervals(Center()),
           crs=SMAPCRS,
           mappedcrs=EPSG(4326),
           dim=X(),
        )
        latlookup = Mapped(latvec;
            order=ReverseOrdered(),
            span=Irregular(latbounds),
            sampling=Intervals(Center()),
            crs=SMAPCRS,
            mappedcrs=EPSG(4326),
            dim=Y(),
        )
        return X(lonlookup), Y(latlookup)
    else
        error("projection $proj not supported")
    end
end

# TODO actually add metadata to the dict
DD.metadata(wrapper::SMAPhdf5) = RA._metadatadict(SMAPsource)

function DD.layerdims(ds::SMAPhdf5)
    keys = RA.cleankeys(RA.layerkeys(ds))
    # All dims are the same
    NamedTuple{keys}(map(_ -> SMAPDIMTYPES, keys))
end

function DD.layermetadata(ds::SMAPhdf5)
    keys = RA.cleankeys(RA.layerkeys(ds))
    NamedTuple{keys}(map(_ -> DD.metadata(ds), keys))
end

Base.keys(ds::SMAPhdf5) = RA.cleankeys(keys(parent(ds)[SMAPGEODATA]))
Base.parent(wrapper::SMAPhdf5) = wrapper.ds
Base.getindex(wrapper::SMAPhdf5, key) = SMAPvar(wrapper.ds[_smappath(key)])

_smappath(key::RA.Key) = SMAPGEODATA * "/" * string(key)

struct SMAPvar{DS} <: AbstractArray{Float32,2}
    ds::DS
end

Base.parent(wrapper::SMAPvar) = wrapper.ds
Base.eltype(wrapper::SMAPvar) = eltype(parent(wrapper))
Base.size(wrapper::SMAPvar) = SMAPSIZE
Base.ndims(wrapper::SMAPvar) = length(SMAPSIZE)
Base.getindex(wrapper::SMAPvar, args...) = getindex(parent(wrapper), args...)
Base.setindex!(wrapper::SMAPvar, args...) = setindex!(parent(wrapper), args...)
Base.Array(wrapper::SMAPvar) = Array(parent(wrapper))
Base.collect(wrapper::SMAPvar) = collect(parent(wrapper))

DA.eachchunk(var::SMAPvar) = DA.GridChunks(var, size(var))
DA.haschunks(var::SMAPvar) = DA.Unchunked()

# Raster ######################################################################

function RA.FileArray(ds::SMAPhdf5, filename::AbstractString; key, kw...)
    RA.FileArray(ds[key], filename; key, kw...)
end
function RA.FileArray(var::SMAPvar, filename::AbstractString; key, kw...)
    T = eltype(var)
    N = ndims(var)
    eachchunk = DA.eachchunk(var)
    haschunks = DA.haschunks(var)
    RA.FileArray{SMAPsource,T,N}(filename, SMAPSIZE; key, eachchunk, haschunks, kw...)
end

function Base.open(f::Function, A::RA.FileArray{SMAPsource}; kw...)
    RA._open(SMAPsource, RA.filename(A); key=RA.key(A), kw...) do var
        f(RA.RasterDiskArray{SMAPsource}(var)) 
    end
end
    
DA.writeblock!(A::RA.RasterDiskArray{SMAPsource}, v, r::AbstractUnitRange...) = A[r...] = v

RA.haslayers(::Type{SMAPsource}) = true

# Stack ########################################################################

function RA.FileStack{SMAPsource}(ds::SMAPhdf5, filename::AbstractString; write=false, keys)
    keys = map(Symbol, keys isa Nothing ? RA.layerkeys(ds) : keys) |> Tuple
    type_size_ec_hc = map(keys) do key
        var = RA.RasterDiskArray{SMAPsource}(ds[key])
        eltype(var), size(var), DA.eachchunk(var), DA.haschunks(var)
    end
    layertypes = map(x->x[1], type_size_ec_hc)
    layersizes = map(x->x[2], type_size_ec_hc)
    eachchunk = map(x->x[3], type_size_ec_hc)
    haschunks = map(x->x[4], type_size_ec_hc)
    RA.FileStack{SMAPsource,keys}(filename, layertypes, layersizes, eachchunk, haschunks, write)
end
function RA.OpenStack(fs::RA.FileStack{SMAPsource,K}; kw...) where K
    ds = HDF5.h5open(RA.filename(fs); kw...)
    RA.OpenStack{SMAPsource,K}(SMAPhdf5(ds))
end
Base.close(os::RA.OpenStack{SMAPsource}) = nothing # HDF5 handles this apparently?

# Series #######################################################################

"""
    smapseries(filenames::AbstractString; kw...)
    smapseries(filenames::Vector{<:AbstractString}, dims=nothing; kw...)

[`RasterSeries`](@ref) loader for SMAP files and whole folders of files,
organised along the time dimension. Returns a [`RasterSeries`](@ref).

# Arguments

- `filenames`: A `String` path to a directory of SMAP files,
    or a vector of `String` paths to specific files.
- `dims`: `Tuple` containing `Ti` dimension for the series.
    Automatically generated form `filenames` unless passed in.

# Keywords

- `kw`: Passed to `RasterSeries`.
"""
function RA.smapseries(dir::AbstractString; kw...)
    RA.smapseries(joinpath.(dir, filter_ext(dir, ".h5")); kw...)
end
function RA.smapseries(filenames::Vector{<:AbstractString}, dims=nothing; kw...)
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
        dims = (_smap_timedim(timeseries),)
    else
        usedpaths = filenames
    end
    # Show errors after all files load, or you can't see them:
    if length(errors) > 0
        println("Some errors thrown during file load: ")
        println.(errors)
    end
    RasterSeries(usedpaths, dims; child=RasterStack, duplicate_first=true, kw...)
end


# Utils ########################################################################

function RA._open(f, ::Type{SMAPsource}, filename::AbstractString; key=nothing, kw...)
    isfile(filename) || RA._filenotfound_error(filename)
    HDF5.h5open(filename; kw...) do ds
        RA._open(f, SMAPsource, SMAPhdf5(ds); key, kw...)
    end
end
function RA._open(f, ::Type{SMAPsource}, ds::SMAPhdf5; key=nothing, kw...)
    RA.cleanreturn(f(key isa Nothing ? ds : ds[key]))
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
function _smap_timedim(times::AbstractVector)
    Ti(Sampled(times, ForwardOrdered(), Regular(Hour(3)), Intervals(Start()), RA._metadatadict(SMAPsource)))
end
