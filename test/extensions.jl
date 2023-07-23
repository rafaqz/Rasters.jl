using Rasters

# These methods should throw a BackendException
@test_throws Rasters.BackendException resample(1, 2, 3)
@test_throws Rasters.BackendException warp(1, 2, 3)

import ArchGDAL

# Should throw a MethodError now that ArchGDAL is loaded
@test_throws MethodError resample(1, 2, 3)
@test_throws MethodError warp(1, 2, 3)