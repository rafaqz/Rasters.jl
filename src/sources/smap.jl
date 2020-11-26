using .HDF5

export SMAPstack, SMAPseries, SMAPdimMetadata, SMAParrayMetadata, SMAPstackMetadata

const SMAPMISSING = -9999.0
const SMAPGEODATA = "Geophysical_Data"
const SMAPCRS = ProjString("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")

"""
    SMAPdimMetadata(val::Union{Dict,NamedTuple})
    SMAPdimMetadata(pairs::Pair...) => SMAPdimMetadata{Dict}
    SMAPdimMetadata(; kw...) => SMAPdimMetadata{NamedTuple}

`DimMetadata` wrapper for `SMAParray` dimensions.
"""
struct SMAPdimMetadata{T} <: AbstractDimMetadata{T}
    val::T
end

"""
    SMAParrayMetadata(val::Union{Dict,NamedTuple})
    SMAParrayMetadata(pairs::Pair...) => SMAParrayMetadata{Dict}
    SMAParrayMetadata(; kw...) => SMAParrayMetadata{NamedTuple}

`ArrayMetadata` wrapper for `SMAParray`.
"""
struct SMAParrayMetadata{T} <: AbstractArrayMetadata{T}
    val::T
end

"""
    SMAPstackMetadata(val::Union{Dict,NamedTuple})
    SMAPstackMetadata(pairs::Pair...) => SMAPstackMetadata{Dict}
    SMAPstackMetadata(; kw...) => SMAPstackMetadata{NamedTuple}

`StackMetadata` wrapper for `SMAPstack`.
"""
struct SMAPstackMetadata{T} <: AbstractStackMetadata{T}
    val::T
end

# Stack ########################################################################

"""
    SMAPstack(filename::String; 
              dims=nothing,
              refdims=nothing,
              window=())

`AbstractGeoStack` for [SMAP](https://smap.jpl.nasa.gov/) datasets.

The simplicity of the format means `dims` and `metadata` are the same for all stack layers,
so we store them as stack fields.

## Arguments

- `filename`: `String` path to a SMAP .h5 file.

## Keyword arguments

- `dims`: Dimensions held on the stack as all layers have identical `Dimension`s. 
  These are loaded from the HDF5 by default, but can be passed in to improve performance,
  as is done by [`SMAPseries`](@ref), 
- `refdims`: As for `dims`. Often the position time `Dimension` from the `SMAPseries`.
- `metadata`: [`SMAPstackMetadata`](@ref) object. As for `dims`.
- `window`: Like `view` but lazy, for disk based data. Can be a tuple of Dimensions,
  selectors or regular indices. These will be applied when the data is loaded or indexed into.
"""
struct SMAPstack{T,D,R,W,M} <: DiskGeoStack{T}
    filename::T
    dims::D
    refdims::R
    window::W
    metadata::M
end
SMAPstack(filename::String;
          dims=smapread(smapdims, filename),
          refdims=(smap_timedim(smap_timefrompath(filename)),),
          window=(),
          metadata=smapread(smapmetadata, filename),
         ) =
    SMAPstack(filename, dims, refdims, window, metadata)


# Dummy SMAParray type. Actually generates a GeoArray
struct SMAParray end

# SMAP has fixed dims for all layers, so we store them on the stack.
dims(stack::SMAPstack, key::Key...) = stack.dims
dims(stack::SMAPstack, dim) = dims(dims(stack), dim)
refdims(stack::SMAPstack) = stack.refdims
metadata(stack::SMAPstack) = stack.metadata
missingval(stack::SMAPstack, key::Key...) = SMAPMISSING
childtype(stack::SMAPstack) = SMAParray
childkwargs(stack::SMAPstack) = ()

withsource(f, ::Type{SMAParray}, path::AbstractString, key...) =
    smapread(f, path)
withsourcedata(f, ::Type{SMAParray}, path::AbstractString, key) =
    smapread(path) do d
        f(d[smappath(string(key))])
    end



# Base methods

# Override getindex as we already have `dims` - they are
# fixed for the whole stack
Base.getindex(s::SMAPstack, key::Key, i1::Integer, I::Integer...) =
    smapread(filename(s)) do file
        dataset = file[smappath(key)]
        _window = maybewindow2indices(dataset, _dims, window(s))
        readwindowed(dataset, _window, I...)
    end
Base.getindex(s::SMAPstack, key::Key, I::Union{Colon,Integer,AbstractArray}...) =
    smapread(filename(s)) do file
        dataset = file[smappath(key)]
        _dims = dims(s)
        _window = maybewindow2indices(dataset, _dims, window(s))
        _dims, _refdims = slicedims(slicedims(_dims, refdims(s), _window)..., I)
        A = readwindowed(dataset, _window, I...)
        GeoArray(A, _dims, _refdims, Symbol(key), metadata(s), missingval(s))
    end

# HDF5 uses `names` instead of `keys` so we have to special-case it
Base.keys(stack::SMAPstack) =
    smapread(filename(stack)) do dataset
        cleankeys(names(dataset[SMAPGEODATA]))
    end


# Series #######################################################################

"""
    SMAPseries(filenames::AbstractString; kwargs...)
    SMAPseries(filenames::Vector{<:AbstractString}, dims=nothing; kwargs...)

[`GeoSeries`](@ref) loader for SMAP files and whole folders of files, 
organised along the time dimension. Returns a [`GeoSeries`](@ref).

## Arguments

- `filenames`: A `String` path to a directory of SMAP files,
  or a vector of `String` paths to specific files.
- `dims`: `Tuple` containing `Ti` dimension for the series. 
  Automatically generated form `filenames` unless passed in.

## Keyword Arguments

- `kwargs`: Passed to `GeoSeries`.
"""
SMAPseries(dir::AbstractString; kwargs...) =
    SMAPseries(joinpath.(dir, filter_ext(dir, ".h5")); kwargs...)
SMAPseries(filenames::Vector{<:AbstractString}, dims=nothing; kwargs...) = begin
    if dims isa Nothing
        usedpaths = String[]
        timeseries = []
        errors = []
        for filename in filenames
            println("Loading SMAP file: $filename")
            try
                t = smap_timefrompath(filename)
                push!(timeseries, t)
                push!(usedpaths, filename)
            catch e
                push!(errors, e)
            end
        end
        # Use the first files time dim as a template, but join vals into an array of times.
        timedim = smap_timedim(timeseries)
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
        dims=smapread(smapdims, first(filenames)),
        metadata=smapread(smapmetadata, first(filenames)),
    )
    GeoSeries(usedpaths, (timedim,); childtype=SMAPstack, childkwargs=childkwargs, kwargs...)
end

Base.:*(hrs::Int, ::Type{T}) where T<:Period = T(hrs)

# Utils ########################################################################

smapread(f, filepath::AbstractString) = h5open(f, filepath)

readwindowed(A::HDF5Dataset, window::Tuple{}) = HDF5.read(A)

smappath(key::Key) = SMAPGEODATA * "/" * string(key)

smap_timefrompath(path::String) = begin
    dateformat = DateFormat("yyyymmddTHHMMSS")
    dateregex = r"SMAP_L4_SM_gph_(\d+T\d+)_"
    datematch = match(dateregex, path)
    if !(datematch === nothing)
        DateTime(datematch.captures[1], dateformat)
    else
        error("Date/time not correctly formatted in path: $path")
    end
end

smap_timedim(t::DateTime) = smap_timedim(t:Hour(3):t)
smap_timedim(times::AbstractVector) =
    Ti(times, mode=Sampled(Ordered(), Regular(Hour(3)), Intervals(Start())))

# TODO actually add metadata to the dict
smapmetadata(dataset::HDF5.HDF5File) = SMAPstackMetadata(Dict())

smapdims(dataset::HDF5.HDF5File) = begin
    proj = read(attrs(root(dataset)["EASE2_global_projection"]), "grid_mapping_name")
    if proj == "lambert_cylindrical_equal_area"
        # There are matrices for lookup but all rows/colums are identical.
        # For performance and simplicity we just take a vector slice for each dim.
        extent = attrs(root(dataset)["Metadata/Extent"])
        lonbounds = read(extent["westBoundLongitude"]), read(extent["eastBoundLongitude"])
        latbounds = read(extent["southBoundLatitude"]), read(extent["northBoundLatitude"])
        latvec = read(root(dataset)["cell_lat"])[1, :]
        lonvec = read(root(dataset)["cell_lon"])[:, 1]
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
        (Lon(lonvec; mode=lonmode), Lat(latvec; mode=latmode))
    else
        error("projection $proj not supported")
    end
end
