using HDF5

export SMAPseries
export smapseries

const SMAPMISSING = -9999.0
const SMAPGEODATA = "Geophysical_Data"
const SMAPCRS = ProjString("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
const SMAPSIZE = (3856, 1624)
const SMAPDIMTYPES = (X, Y) 

# Dataset wrapper ###############################################################
# Becauase SMAP is just one of manyu HDF5 formats,
# we wrap it in SMAPhdr5 and SMAPvar wrappers

struct SMAPhdf5{T}
    ds::T
end
filename(wrapper::SMAPhdf5) = wrapper.filename
Base.parent(wrapper::SMAPhdf5) = wrapper.ds
Base.getindex(wrapper::SMAPhdf5, path) = wrapper.ds[path]

_smappath(key::Key) = SMAPGEODATA * "/" * string(key)

struct SMAPvar{T}
    ds::T
end
Base.parent(wrapper::SMAPvar) = wrapper.ds

# GeoArray ######################################################################

function FileArray(var::SMAPvar, filename::AbstractString; kw...)
    T = eltype(parent(var))
    N = length(SMAPSIZE)
    FileArray{_SMAP,T,N}(filename, SMAPSIZE; kw...)
end

function Base.open(f::Function, A::FileArray{_SMAP})
    _read(ds -> f(ds), _SMAP, filename(A); key=key(A))
end

# Stack ########################################################################

hasstackfile(::Type{_SMAP}) = true

function Base.getindex(fs::FileStack{_SMAP}, key)
   _read(_SMAP, filename(fs); key) do var
       FileArray(SMAPvar(var), filename(fs); key)
   end
end

Base.keys(ds::SMAPhdf5) = cleankeys(keys(parent(ds)[SMAPGEODATA]))

@inline function DD.layerdims(ds::SMAPhdf5)
    keys = cleankeys(layerkeys(ds))
    # All dims are the same
    NamedTuple{keys}(map(_ -> SMAPDIMTYPES, keys))
end

@inline function DD.layermetadata(ds::SMAPhdf5)
    keys = cleankeys(layerkeys(ds))
    NamedTuple{keys}(map(_ -> DD.metadata(ds), keys))
end

# TODO actually add metadata to the dict
DD.metadata(wrapper::SMAPhdf5) = Metadata{_SMAP}(Dict())

missingval(ds::SMAPhdf5) = SMAPMISSING
layermissingval(ds::SMAPhdf5) = SMAPMISSING
layerkeys(ds::SMAPhdf5) = keys(ds)
layersizes(ds::SMAPhdf5, keys) = map(_ -> SMAPSIZE, keys)

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
    SMAPseries(joinpath.(dir, filter_ext(dir, ".h5")); kw...)
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
    dims, metadata =_read(_SMAP, first(filenames)) do ds
        DD.dims(ds), DD.metadata(ds)
    end
    childkwargs = (; dims, metadata)
    GeoSeries(usedpaths, (timedim,); child=stack, childkwargs=childkwargs, kw...)
end

@deprecate SMAPseries(args...; kw...) smapseries(args...; kw...)

function DD.dims(wrapper::SMAPhdf5)
    dataset = parent(wrapper)
    proj = read(HDF5.attributes(HDF5.root(dataset)["EASE2_global_projection"]), "grid_mapping_name")
    if proj == "lambert_cylindrical_equal_area"
        # There are matrices for lookup but all rows/colums are identical.
        # For performance and simplicity we just take a vector slice for each dim.
        extent = HDF5.attributes(HDF5.root(dataset)["Metadata/Extent"])
        lonbounds = read(extent["westBoundLongitude"]), read(extent["eastBoundLongitude"])
        latbounds = read(extent["southBoundLatitude"]), read(extent["northBoundLatitude"])
        latvec = read(HDF5.root(dataset)["cell_lat"])[1, :]
        lonvec = read(HDF5.root(dataset)["cell_lon"])[:, 1]
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

DD.refdims(wrapper::SMAPhdf5, filename) = (_smap_timedim(_smap_timefromfilename(filename)),)

# Utils ########################################################################

function _read(f, ::Type{_SMAP}, filepath::AbstractString; key=nothing, kw...)
    if key isa Nothing
        h5open(ds -> f(SMAPhdf5(ds)), filepath; kw...)
    else
        h5open(ds -> f(SMAPhdf5(ds)[_smappath(key)]), filepath)
    end
end

function _smap_timefromfilename(filename::String)
    dateformat = DateFormat("yyyymmddTHHMMSS")
    dateregex = r"SMAP_L4_SM_gph_(\d+T\d+)_"
    datematch = match(dateregex, filename)
    if !(datematch === nothing)
        DateTime(datematch.captures[1], dateformat)
    else
        error("Date/time not correctly formatted in path: $filenampathe")
    end
end

_smap_timedim(t::DateTime) = _smap_timedim(t:Hour(3):t)
_smap_timedim(times::AbstractVector) =
    Ti(times, mode=Sampled(Ordered(), Regular(Hour(3)), Intervals(Start())))
