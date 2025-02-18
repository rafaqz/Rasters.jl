# # Writing rasters to Disk

# Rasters.jl supports writing datasets to disk in several different formats. `.grd` files are 
# supported with no extra steps, but all others require loading an extension.
# All other formats are supported with the ArchGDAL and NCDatasets extensions.
# Datasets can (with some caveats) be written to different formats than they were loaded 
# in as, providing file-type conversion for spatial data.

# Some metadata may be lost in formats that store little metadata, or where
# metadata conversion has not been completely implemented.

# # Create a Raster to save
using Rasters

# specify the coodinate reference system [crs]. 
proj = ProjString("+proj=longlat +datum=WGS84 +type=crs")
# or EPSG(4326) or some ESRIWellKnownText string - see GeoFormatTypes.jl for more info!

# define coodinate dimensions
lon, lat = X(25:1:30), Y(25:1:30)

# create a Raster 
ras = Raster(rand(lon, lat); crs=proj) # this generates random numbers with the dimensions given

# # Write Raster to a `.grd` file
# Note: that even though grd is a build in format that should not imply that it is the
# most suitable format for saving your data.
filename = tempname() * ".grd"
write(filename, ras)

# # Write Raster using GDAL (requires ArchGDAL.jl)

# write(filename::AbstractString, A::AbstractRaster; kw...)
#   Write an AbstractRaster to file, guessing the backend from the file extension.
#   Keyword arguments are passed to the write method for the backend.

using Pkg;
Pkg.add("ArchGDAL");
using ArchGDAL

# write Raster as a GeoTIFF
filename = tempname() * ".tif"
write(filename, ras)

# write Raster as a cog (cloud optimized geotiff) by specifing the COG `driver``
filename = tempname() * ".tif"
write(filename, ras, driver="COG")

# you can also pass options to GDAL by passing a dictionary of options. See GDAL 
# specific driver options. [See options specific to the GOG diver.](https://gdal.org/drivers/raster/cog.html)
# `force=true` is required to overwirte an exisiting file
write(filename, ras, driver="COG", options=Dict("BLOCKSIZE" => "512", "RESAMPLING" => "NEAREST"), force=true)

# See the [full list of GDAL drivers](https://gdal.org/drivers/raster/index.html) 
# for other supported formats. 

# # write Raster to a NetCDF file [requires adding NCDatasets.jl]
Pkg.add("NCDatasets");
using NCDatasets

filename = tempname() * ".nc"
write(filename, ras)

# # write Raster to a JDL2 file. 

# JLD2 saves and loads Julia data structures in a format comprising a subset of HDF5, 
# without any dependency on the HDF5 C library. The save and load functions, provided 
# by FileIO, provide a mechanism to read and write data from a JLD2 file.
Pkg.add("FileIO");
using FileIO

filename = tempname() * ".jld2"
save(filename, ras)

