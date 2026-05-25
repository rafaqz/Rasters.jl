using Test, Rasters

# These methods should throw a BackendException
@test_throws Rasters.BackendException resample(1, 2, 3)
@test_throws Rasters.BackendException warp(1, 2, 3)

@test_throws "import ZarrDatasets" Raster("notafile.zarr")
@test_throws "import ArchGDAL" Raster("notafile.tif")
@test_throws "import NCDatasets" Raster("notafile.nc")
@test_throws "import GRIBDatasets" Raster("notafile.grib")

import ArchGDAL

# Should throw a MethodError now that ArchGDAL is loaded
@test_throws MethodError resample(1, 2, 3)
@test_throws MethodError warp(1, 2, 3)
