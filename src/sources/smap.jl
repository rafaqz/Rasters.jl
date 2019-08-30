# A data source for SMAP (Soil moisture active-passive) datatsets.

using HDF5, Dates

export SMAPstack, SMAPseries

const SMAPMISSING = -9999.0
const SMAPEXTENT = "Metadata/Extent"
const SMAPGEODATA = "Geophysical_Data"

struct SMAPstack{T<:AbstractString} <: AbstractGeoStack{T}
    data::T
end

dims(stack::SMAPstack, ds, key::Key) = smapcoords(ds)
refdims(stack::SMAPstack) = (run(smaptime, stack, source(stack)),)
metadata(stack::SMAPstack, args...) = nothing
missingval(stack::SMAPstack, args...) = SMAPMISSING

@inline run(f, stack::SMAPstack, path::AbstractString) = h5open(f, path)

data(s::SMAPstack, ds, key::Key) = 
    GeoArray(read(ds[smappath(key)]), dims(s, key), refdims(s), metadata(s), missingval(s))
data(s::SMAPstack, ds, key::Key, I...) = begin
    data = ds[smappath(key)][I...]
    GeoArray(data, slicedims(dims(s, key), refdims(s), I)..., metadata(s), missingval(s))
end
data(s::SMAPstack, ds, key::Key, I::Vararg{Integer}) = ds[smappath(key)][I...]

# HDF5 uses `names` instead of `keys` so we have to special case it
Base.keys(stack::SMAPstack) = run(ds -> names(ds[SMAPGEODATA]), stack, source(stack))
Base.copy!(dst::AbstractArray, src::SMAPstack, key) = 
    run(ds -> dst .= read(ds[smappath(key)]), src, source(src))

"""
Series loader for SMAP folders (files in the time dimension). 
It outputs a GeoSeries
"""
SMAPseries(path::AbstractString) = SMAPseries(joinpath.(path, filter_ext(path, ".h5")))
SMAPseries(filepaths::Vector{<:AbstractString}, dims=smapseriestime(filepaths); 
           refdims=(), metadata=nothing, childtype=SMAPstack) = 
    GeoSeries(filepaths, dims, refdims, metadata, childtype)

# smap utilities

smappath(key) = joinpath(GeoData.SMAPGEODATA, string(key))

smaptime(ds) = begin
    meta = attrs(root(ds)["time"])
    units = read(meta["units"])
    datestart = replace(replace(units, "seconds since " => ""), " " => "T")
    dt = DateTime(datestart) + Dates.Second(read(meta["actual_range"])[1])
    Time(dt)
end

smapseriestime(filepaths) = begin
    timeseries = h5open.(ds -> smaptime(ds), filepaths)
    timemeta = metadata(first(timeseries))
    (Time(val.(timeseries), timemeta),)
end

smapcoords(ds) = begin
    proj = read(attrs(root(ds)["EASE2_global_projection"]), "grid_mapping_name")
    if proj == "lambert_cylindrical_equal_area"
        # There are matrices for lookup but all rows/colums are identical.
        # For performance we just take a vector slice of each dim.
        latvec = read(root(ds)["cell_lat"])[1, :]
        lonvec = read(root(ds)["cell_lon"])[:, 1]
        (Lon(lonvec), Lat(latvec))
    else
        error("projection $proj not supported")
    end
end
