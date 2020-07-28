
using GeoData, NCDatasets, Statistics, Plots, ArchGDAL

using GeoData, NCDatasets, Plots
pyplot()
filename = download("https://www.unidata.ucar.edu/software/netcdf/examples/tos_O1_2001-2002.nc", "tos_O1_2001-2002.nc")

stack = NCDstack(filename)
keys(stack) # (:tos,)
A = stack[:tos]

A[Ti(1:3:12)] |> plot
savefig("tos_4.png")

A[Ti(Contains(DateTime360Day(2001, 01, 17))), Lat(Between(0, -50)), Lon(Between(100, 160))] |> plot
savefig("tos_jan_australia.png")

mean(A; dims=Ti) |> plot
savefig("mean_tos.png")

A[Lat(Contains(20)), Ti(1)] |> plot
savefig("tos_20deg_lat.png")


geturl(url) = begin
    fname = splitdir(url)[2]
    isfile(fname) || download(url, fname)
    fname
end

write("mean_tos.ncd", NCDarray, mean(A; dims=Ti))

# Load some layers from NetCDF #############################################

# Sea surface temperatures collected by PCMDI for use by the IPCC.
ncurl = "https://www.unidata.ucar.edu/software/netcdf/examples/tos_O1_2001-2002.nc"
ncfilename = geturl(ncurl)

# Create a stack from a netcdf file
stack = NCDstack(ncfilename)
A = first(stack)
A[Ti(1)] |> plot

# Create a GeoSeries from multiple files.
# This uses the same data three times to avoid downloads, you 
# would really use a series of datasets matching the time dimension dates.
dimz = (Ti<|[DateTime360Day(2001, 01, 1), DateTime360Day(2001, 02, 1), DateTime360Day(2001, 03, 1)],)
filenames = [ncfilename, ncfilename, ncfilename]
series = GeoSeries(filenames, dimz; childtype=NCDstack)

# Get a single array from the series
a = series[Near<|DateTime360Day(2001, 01, 1)]["tos"][Lon<|Between(50, 200)]

# Do some manipulations and plotting ###################################

# Plot single slices
view(a, Ti(4)) |> plot

# Or reduce over some dimension
mean(a; dims=Ti) |> plot
std(a; dims=Ti) |> plot

# Other things work too
minimum(replace_missing(a, NaN); dims=Ti) |> plot
maximum(replace_missing(a, NaN); dims=Ti) |> plot
reduce(+, a; dims=Ti) |> plot

# Plot the mean sea surface temperature for australia in the second half of 2002 
stack = NCDstack(ncfilename)
t = Ti<|Between(DateTime360Day(2002, 07, 1), DateTime360Day(2002, 012, 30)) 
stack["tos"][t, Lat<|Between(-45, 0.5), Lon<|Between(110, 160)] |> x->mean(x; dims=Ti) |> plot

# Permute the dimensions in the underlying data
# It stil plots the right way up
permutedims(a, (Lat, Lon, Ti))[Ti(1:3:12)] |> plot

# Line plots have (kind of) useful labels
a = series[1]["tos"]
replace(a[Lat(80)], missing=>NaN) |> plot
a[Lon(170), Ti(10)] |> plot
a[Lat(1:80), Lon(170), Ti(10)] |> plot


# Load a tif with GDAL ######################################################
gdal_url = "https://download.osgeo.org/geotiff/samples/gdal_eg/cea.tif"
filepath = geturl(gdal_url)

array = GDALarray(filepath)
array[Band(1)] |> plot

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

# Move it to a regular memory based GeoArray
memarray = GeoArray(diskarray) 
size(memarray)
dims(memarray)
memarray |> plot
