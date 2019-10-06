using GeoData, Test, Statistics, Dates
using NCDatasets, HDF5, ArchGDAL
using GeoData: Time, formatdims, data, dims2indices, rebuild, window, name

# Loader for external sources
geturl(url) = begin
    fname = splitdir(url)[2]
    isfile(fname) || download(url, fname)
    fname
end

include("array.jl")
include("stack.jl")
include("series.jl")
include("gdal.jl")
include("ncdatasets.jl")
include("smap.jl")
