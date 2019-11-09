using .HDF5

export SMAPstack, SMAPseries

const SMAPMISSING = -9999.0
const SMAPEXTENT = "Metadata/Extent"
const SMAPGEODATA = "Geophysical_Data"


# Stack ########################################################################

@GeoStackMixin struct SMAPstack{} <: AbstractGeoStack{T} end

SMAPstack(data::String; dims=smapapply(smapdims, data), 
          refdims=(), window=(), metadata=ncapply(metadata, data)) =
    SMAPstack(data, dims, refdims, window, metadata)

dims(stack::SMAPstack) = stack.dims
dims(stack::SMAPstack, ::Key) = stack.dims
refdims(stack::SMAPstack) = (safeapply(smaptime, stack, source(stack)),)
metadata(stack::SMAPstack, args...) = nothing
missingval(stack::SMAPstack, args...) = SMAPMISSING

@inline safeapply(f, ::SMAPstack, path::AbstractString) = smapapply(f, path)

data(s::SMAPstack, dataset, key::Key) =
    GeoArray(read(dataset[smappath(key)]), dims(s), refdims(s), metadata(s), missingval(s), Symbol(key))
data(s::SMAPstack, dataset, key::Key, I...) = begin
    data = dataset[smappath(key)][I...]
    GeoArray(data, slicedims(dims(s), refdims(s), I)..., metadata(s), missingval(s), Symbol(key))
end
data(::SMAPstack, dataset, key::Key, I::Vararg{Integer}) = dataset[smappath(key)][I...]

# HDF5 uses `names` instead of `keys` so we have to special-case it
Base.keys(stack::SMAPstack) = 
    safeapply(dataset -> Tuple(Symbol.(names(dataset[SMAPGEODATA]))), stack, source(stack))
Base.copy!(dst::AbstractArray, src::SMAPstack, key) =
    safeapply(dataset -> copy!(dst, dataset[smappath(key)][window2indices(src)...]), src, source(src))
Base.copy!(dst::AbstractGeoArray, src::SMAPstack, key) =
    safeapply(dataset -> copy!(parent(dst), dataset[smappath(key)][window2indices(src)...]), src, source(src))

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

smapapply(f, path) = h5open(f, path)

smappath(key) = joinpath(GeoData.SMAPGEODATA, string(key))

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

smapdims(dataset) = begin
    proj = read(attrs(root(dataset)["EASE2_global_projection"]), "grid_mapping_name")
    if proj == "lambert_cylindrical_equal_area"
        # There are matrices for lookup but all rows/colums are identical.
        # For performance and simplicity we just take a vector slice for each dim.
        latvec = read(root(dataset)["cell_lat"])[1, :]
        lonvec = read(root(dataset)["cell_lon"])[:, 1]
        latgrid=AllignedGrid(order=Ordered(Reverse(), Reverse()))
        (Lon(lonvec), Lat(latvec; grid=latgrid))
    else
        error("projection $proj not supported")
    end
end
