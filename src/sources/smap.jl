using HDF5

export SMAP, SMAPstack, SMAPseries

const SMAPMISSING = -9999.0
const SMAPGEODATA = "Geophysical_Data"
const SMAPCRS = ProjString("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

hasstackfile(::Type{_SMAP}) = true

struct SMAPhdf5{T}
    ds::T
end
Base.parent(wrapper::SMAPhdf5) = wrapper.ds

# Stack ########################################################################

"""
    SMAPstack <: DiskGeoStack

    SMAPstack(filename::String; dims=nothing, refdims=nothing, window=())

`AbstractGeoStack` for [SMAP](https://smap.jpl.nasa.gov/) datasets.

The simplicity of the format means `dims` and `metadata` are the same for all stack layers,
so we store them as stack fields.

# Arguments

- `filename`: `String` path to a SMAP .h5 file.

# Keywords

- `dims`: Dimensions held on the stack as all layers have identical `Dimension`s.
    These are loaded from the HDF5 by default, but can be passed in to improve performance,
    as is done by [`SMAPseries`](@ref),
- `refdims`: As for `dims`. Often the position time `Dimension` from the `SMAPseries`.
- `metadata`: `Metadata` object. As for `dims`.
- `window`: Like `view` but lazy, for disk based data. Can be a tuple of Dimensions,
    selectors or regular indices. These will be applied when the data is loaded or indexed.
"""
struct SMAPstack end



function Base.getindex(fs::FileStack{_SMAP}, key)
   _read(_SMAP, filename(fs)) do ds
       var = ds[string(key)]
       FileArray(var, key, filename(fs))
   end
end


# SMAP has fixed dims for all layers, so we store them on the stack.
DD.dims(stack::SMAPstack, dim) = dims(dims(stack), dim)
DD.dims(stack::SMAPstack, key::Key...) = stack.dims
DD.refdims(stack::SMAPstack) = stack.refdims
DD.metadata(stack::SMAPstack) = stack.metadata
missingval(stack::SMAPstack, key::Key...) = SMAPMISSING
childtype(stack::SMAPstack) = SMAParray
childkwargs(stack::SMAPstack) = ()

# Dummy SMAParray type. Actually generates a GeoArray
struct SMAParray end

withsource(f, ::Type{SMAParray}, path::AbstractString, key...) = _smapread(f, path)
withsourcedata(f, ::Type{SMAParray}, path::AbstractString, key) =
    _smapread(path) do d
        f(d[_smappath(string(key))])
    end


# Base methods

# Override getindex as we already have `dims` - they are
# fixed for the whole stack
function Base.getindex(s::FileStack{_SMAP}, key::Key, i1::Integer, I::Integer...)
    _smapread(filename(s)) do wrapper
        dataset = parent(wrapper)
        smapdata = file[_smappath(key)]
        _window = maybewindow2indices(smapdata, _dims, window(s))
        readwindowed(smapdata, _window, I...)
    end
end
function Base.getindex(s::SMAPstack, key::Key)
    _smapread(filename(s)) do wrapper
        dataset = parent(wrapper)
        smapdata = dataset[_smappath(key)]
        dims_ = dims(s)
        window_ = maybewindow2indices(smapdata, dims_, window(s))
        dims_, refdims_ = DD.slicedims(DD.slicedims(dims_, refdims(s), window_)..., I)
        A = FileArray(smapdata, )
        GeoArray(A, dims_, refdims_, Symbol(key), metadata(s), missingval(s))
    end
end

# HDF5 uses `names` instead of `keys` so we have to special-case it
function Base.keys(stack::SMAPstack)
    _smapread(filename(stack)) do dataset
        cleankeys(keys(dataset[SMAPGEODATA]))
    end
end

layerkeys(stack::SMAPstack) = keys(stack)


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
function SMAPseries(dir::AbstractString; kw...)
    SMAPseries(joinpath.(dir, filter_ext(dir, ".h5")); kw...)
end
function SMAPseries(filenames::Vector{<:AbstractString}, dims=nothing; kw...)
    if dims isa Nothing
        usedpaths = String[]
        timeseries = []
        errors = []
        for filename in filenames
            try
                t = _smap_timefrompath(filename)
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
    childkwargs = (
        dims=_smapread(_smapdims, first(filenames)),
        metadata=_smapread(_smapmetadata, first(filenames)),
    )
    GeoSeries(usedpaths, (timedim,); childtype=SMAPstack, childkwargs=childkwargs, kw...)
end

Base.:*(hrs::Int, ::Type{T}) where T<:Period = T(hrs)

readwindowed(A::HDF5.Dataset, window::Tuple{}) = HDF5.read(A)


# Utils ########################################################################

_smapread(f, filepath::AbstractString) = _read(f, _SMAP, filepath)

_read(f, ::Type{_SMAP}, filepath::AbstractString) = h5open(ds -> f(SMAPhdf5(ds)), filepath)

_smappath(key::Key) = SMAPGEODATA * "/" * string(key)

function _smap_timefrompath(path::String)
    dateformat = DateFormat("yyyymmddTHHMMSS")
    dateregex = r"SMAP_L4_SM_gph_(\d+T\d+)_"
    datematch = match(dateregex, path)
    if !(datematch === nothing)
        DateTime(datematch.captures[1], dateformat)
    else
        error("Date/time not correctly formatted in path: $path")
    end
end

_smap_timedim(t::DateTime) = _smap_timedim(t:Hour(3):t)
_smap_timedim(times::AbstractVector) =
    Ti(times, mode=Sampled(Ordered(), Regular(Hour(3)), Intervals(Start())))

# TODO actually add metadata to the dict
DD.metadata(wrapper::HDF5.File) = Metadata{:SMAP}(Dict())

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
