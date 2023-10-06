# # Writing rasters to Disk

# Rasters.jl supports writing datasets to disk in several differnt formats. grd files are 
# currently the only supported file format that does not require loading an extension. 
# All other formats are supported with the ArchGDAL, HDF5, and NCDatasets extensions.
# Datasets can (with some caveats) be written to different formats than they were loaded 
# in as, providing file-type conversion for spatial data.

# Some metadata may be lost in formats that store little metadata, or where
# metadata conversion has not been completely implemented.

# # create a Raster dataset for saving
using Rasters

# specify the coodinate reference system [crs]. 
proj = "+proj=longlat +datum=WGS84 +no_defs +type=crs"

# define coodinate dimensions
lon, lat = X(25:1:30), Y(25:1:30)

# create a Raster 
ras = Raster(rand(lon, lat); crs=proj) # this generates random numbers with the dimensions given

# # write Raster to a grd file
# Note: that even though grd is a build in format that should not imply that it is the
# most suitable format for saving your data.
filename = tempname() * ".grd"
write(filename, ras)

# # write Raster using GDAL[requires adding ArchGDAL.jl]

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

# # write Raster to a HDF5 file [requires adding HDF5.jl]
using Pkg;
Pkg.add("HDF5");
using HDF5

filename = tempname() * ".hdf"
write(filename, ras)