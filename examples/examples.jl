# First time:
# ] dev http://github.com/rafaqz/CoordinateReferenceSystemsBase.jl
# ] dev http://github.com/rafaqz/DimensionalData.jl
# ] dev http://github.com/rafaqz/GeoData.jl
# Load data #######################################################################

using GeoData, NCDatasets, Statistics, Plots, ArchGDAL
using GeoData: Time
geturl(url) = begin
    fname = splitdir(url)[2]
    isfile(fname) || download(url, fname)
    fname
end

pyplot()


# Load some layers from NetCDF #############################################

# Sea surface temperatures collected by PCMDI for use by the IPCC.
ncurl = "https://www.unidata.ucar.edu/software/netcdf/examples/tos_O1_2001-2002.nc"
ncfilename = geturl(ncurl)

# Create a stack from a netcdf file
stack = NCstack(ncfilename)

# Create a GeoSeries from multiple files.
# This uses the same data three times to avoid downloads, you 
# would really use a series of datasets matching the time dimension dates.
dimz = (Time<|[DateTime360Day(2001, 01, 1), DateTime360Day(2001, 02, 1), DateTime360Day(2001, 03, 1)],)
filenames = [ncfilename, ncfilename, ncfilename]
series = GeoSeries(filenames, dimz; childtype=NCstack)

# Get a single array from the series
a = series[Near<|DateTime360Day(2001, 01, 1)]["tos"][Lon<|Between(50, 200)]
dims(a)


# Do some more manipulations and plotting ###################################

# Plot the mean sea surface temperature for australia in the second half of 2002 
stack = NCstack(ncfilename)
t = Time<|Between(DateTime360Day(2002, 07, 1), DateTime360Day(2002, 012, 30)) 
stack["tos"][t, Lat<|Between(-45, 0.5), Lon<|Between(110, 160)] |> x->mean(x; dims=Time) |> plot

# Plot single slices
view(a, Time(4)) |> plot

# Or reduce over some dimension
mean(a; dims=Time) |> plot
std(a; dims=Time) |> plot

# Other things work too
minimum(x->ismissing(x) ? NaN : x, a; dims=Time) |> plot
maximum(x->ismissing(x) ? NaN : x, a; dims=Time) |> plot
reduce(+, a; dims=Time) |> plot

# Permute the dimensions in the underlying data
# It stil plots the right way up
permutedims(a, (Lat, Lon, Time)) |> plot

# Line plots have (kind of) useful labels
a = series[1]["tos"]
replace(a[Lat(80)], missing=>NaN) |> plot
a[Lon(170), Time(10)] |> plot
a[Lat(1:80), Lon(170), Time(10)] |> plot


# Load a tif with GDAL ######################################################
gdal_url = "https://download.osgeo.org/geotiff/samples/gdal_eg/cea.tif"
filepath = geturl(gdal_url)

array = GDALarray(filepath)
GeoData.window(array)
array[1:2, 3:4, 1]
dims(array)

filepaths = [filepath, filepath, filepath]
keyz = (:one, :two, :three)
diskstack = GDALstack(filepaths, keyz)
memstack = GeoStack(diskstack)
# Plot a section of the tif file
diskstack[:two][Band(1), Lat(1:100), Lon(1:150)] |> plot

# Get an array from the stack, 
diskarray = diskstack[:two, Band(1), Lat(1:100), Lon(1:150)]; # show is broken for GDALarray
# The diskarray is subsetted, like SubArray, but without an actual array underneath yet
GeoData.window(diskarray)
size(diskarray)
metadata(diskarray)

# Move it to a regular disk based GeoArray
memarray = GeoArray(diskarray) 
size(memarray)
dims(memarray)
memarray |> plot
