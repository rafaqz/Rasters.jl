using .HDF5

export SMAPstack, SMAPseries

const SMAPMISSING = -9999.0
const SMAPEXTENT = "Metadata/Extent"
const SMAPGEODATA = "Geophysical_Data"


struct SMAPmetadata{K,V} <: AbstractArrayMetadata{K,V}
    val::Dict{K,V}
end

struct SMAPdimMetadata{K,V} <: AbstractDimMetadata{K,V}
    val::Dict{K,V}
end

# Stack ########################################################################

struct SMAPstack{T,D,R,W,M} <: DiskGeoStack{T}
    filename::T
    dims::D
    refdims::R
    window::W
    metadata::M
end

SMAPstack(filename::String; 
          dims=smapapply(smapdims, filename),
          refdims=(smapapply(smaptime, filename),), 
          window=(), 
          metadata=smapapply(smapmetadata, filename)) =
    SMAPstack(filename, dims, refdims, window, metadata)

dims(stack::SMAPstack, ::Key) = stack.dims
refdims(stack::SMAPstack) = stack.refdims
metadata(stack::SMAPstack) = stack.metadata
missingval(stack::SMAPstack) = SMAPMISSING

@inline safeapply(f, ::SMAPstack, path::AbstractString) = smapapply(f, path)

@inline Base.getindex(s::SMAPstack, key::Key, i1::Integer, I::Integer...) =
    smapapply(filename(s)) do file
        dataset = file[smappath(key)]
        _window = maybewindow2indices(dataset, _dims, window(s))
        smapread(dataset, _window, I...)
    end
@inline Base.getindex(s::SMAPstack, key::Key, I::Union{Colon,Integer,AbstractArray}...) =
    smapapply(filename(s)) do file
        dataset = file[smappath(key)]
        _dims = dims(s)
        _window = maybewindow2indices(dataset, _dims, window(s))
        _dims, _refdims = slicedims(slicedims(_dims, refdims(s), _window)..., I)
        A = smapread(dataset, _window, I...)
        GeoArray(A, _dims, _refdims, metadata(s), missingval(s), string(key))
    end

# HDF5 uses `names` instead of `keys` so we have to special-case it
Base.keys(stack::SMAPstack) =
    smapapply(filename(stack)) do dataset 
        Tuple(Symbol.(names(dataset[SMAPGEODATA]))) 
    end

Base.copy!(dst::AbstractArray, src::SMAPstack, key) =
    smapapply(filename(src)) do file
        dataset = file[smappath(key)]
        _window = maybewindow2indices(dataset, _dims, window(dst))
        copy!(dst, smapread(dataset, _window))
    end
Base.copy!(dst::AbstractGeoArray, src::SMAPstack, key) =
    smapapply(filename(src)) do file
        dataset = file[smappath(key)]
        _window = maybewindow2indices(dataset, dims(src), window(src))
        copy!(parent(dst), smapread(dataset, _window))
    end

# Series #######################################################################

"""
Series loader for SMAP folders (files in the time dimension).
It outputs a GeoSeries
"""
SMAPseries(path::AbstractString; kwargs...) =
    SMAPseries(joinpath.(path, filter_ext(path, ".h5")); kwargs...)
SMAPseries(filepaths::Vector{<:AbstractString}, dims=smapseriestime(filepaths);
           childtype=SMAPstack,
           childdims=h5open(smapdims, first(filepaths)),
           window=(), kwargs...) =
    GeoSeries(filepaths, dims; childtype=childtype, childdims=childdims, window=window, kwargs...)


# Utils ########################################################################

smapapply(f, filepath) = h5open(f, filepath)

smapread(s::SMAPstack, key, window, I...) =
    smapapply(filename(s)) do file 
        smapread(file[smappath(key)], window, I...) 
    end
smapread(A, window::Tuple{}) = HDF5.read(A)
smapread(A, window::Tuple{}, I...) = A[I...]
smapread(A, window, I...) = A[Base.reindex(window, I)...]
smapread(A, window) = A[window...]


smappath(key) = joinpath(SMAPGEODATA, string(key))

smaptime(dataset) = begin
    meta = attrs(root(dataset)["time"])
    units = read(meta["units"])
    datestart = replace(replace(units, "seconds since " => ""), " " => "T")
    dt = DateTime(datestart) + Dates.Second(read(meta["actual_range"])[1])
    Time(dt)
end

smapseriestime(filepaths) = begin
    timeseries = h5open.(dataset -> smaptime(dataset), filepaths)
    timemeta = metadata(first(timeseries))
    (Time(val.(timeseries); metadata=timemeta),)
end

smapmetadata(dataset) = SMAPmetadata(Dict())

smapdims(dataset) = begin
    proj = read(attrs(root(dataset)["EASE2_global_projection"]), "grid_mapping_name")
    if proj == "lambert_cylindrical_equal_area"
        # There are matrices for lookup but all rows/colums are identical.
        # For performance and simplicity we just take a vector slice for each dim.
        latvec = read(root(dataset)["cell_lat"])[1, :]
        lonvec = read(root(dataset)["cell_lon"])[:, 1]
        latgrid=AllignedGrid(order=Ordered(Reverse(), Reverse(), Forward()))
        (Lon(lonvec), Lat(latvec; grid=latgrid))
    else
        error("projection $proj not supported")
    end
end
