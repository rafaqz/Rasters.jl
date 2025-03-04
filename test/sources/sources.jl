using Rasters, Test
using Rasters: sourcetrait

@test sourcetrait("x", ".nc") == Rasters.NCDsource()
@test sourcetrait("x", ".nc4") == Rasters.NCDsource()
@test sourcetrait("x", ".h5") == Rasters.NCDsource()
@test sourcetrait("x", ".grd") == Rasters.GRDsource()
@test sourcetrait("x", ".gri") == Rasters.GRDsource()
@test sourcetrait("x", ".grib") == Rasters.GRIBsource()
@test sourcetrait("x", ".tif") == Rasters.GDALsource()
@test sourcetrait("x", ".zarr") == Rasters.Zarrsource()

@test sourcetrait("x", :netcdf) == Rasters.NCDsource()
@test sourcetrait("x", :grd) == Rasters.GRDsource()
@test sourcetrait("x", :grib) == Rasters.GRIBsource()
@test sourcetrait("x", :gdal) == Rasters.GDALsource()
@test sourcetrait("x", :zarr) == Rasters.Zarrsource()

@test sourcetrait("x", Rasters.NCDsource()) == Rasters.NCDsource()
@test sourcetrait("x", Rasters.GRDsource()) == Rasters.GRDsource()
@test sourcetrait("x", Rasters.GRIBsource()) == Rasters.GRIBsource()
@test sourcetrait("x", Rasters.GDALsource()) == Rasters.GDALsource()
@test sourcetrait("x", Rasters.Zarrsource()) == Rasters.Zarrsource()

@test sourcetrait("x", Rasters.NCDsource) == Rasters.NCDsource()
@test sourcetrait("x", Rasters.GRDsource) == Rasters.GRDsource()
@test sourcetrait("x", Rasters.GRIBsource) == Rasters.GRIBsource()
@test sourcetrait("x", Rasters.GDALsource) == Rasters.GDALsource()
@test sourcetrait("x", Rasters.Zarrsource) == Rasters.Zarrsource()

@test sourcetrait("x", :netcdf) == Rasters.NCDsource()
@test sourcetrait("x", :grd) == Rasters.GRDsource()
@test sourcetrait("x", :grib) == Rasters.GRIBsource()
@test sourcetrait("x", :gdal) == Rasters.GDALsource()
@test sourcetrait("x", :zarr) == Rasters.Zarrsource()

@test sourcetrait(".nc") == Rasters.NCDsource()
@test sourcetrait(".nc4") == Rasters.NCDsource()
@test sourcetrait(".h5") == Rasters.NCDsource()
@test sourcetrait(".grd") == Rasters.GRDsource()
@test sourcetrait(".grib") == Rasters.GRIBsource()
@test sourcetrait(".tif") == Rasters.GDALsource()
@test sourcetrait(".zarr") == Rasters.Zarrsource()
@test sourcetrait(".zarr/") == Rasters.Zarrsource()

@test sourcetrait("x.nc") == Rasters.NCDsource()
@test sourcetrait("x.nc4") == Rasters.NCDsource()
@test sourcetrait("x.h5") == Rasters.NCDsource()
@test sourcetrait("x.grd") == Rasters.GRDsource()
@test sourcetrait("x.gri") == Rasters.GRDsource()
@test sourcetrait("x.grib") == Rasters.GRIBsource()
@test sourcetrait("x.tif") == Rasters.GDALsource()
@test sourcetrait("x.zarr") == Rasters.Zarrsource()
@test sourcetrait("x.zarr/") == Rasters.Zarrsource()

@test sourcetrait(:netcdf) == Rasters.NCDsource()
@test sourcetrait(:grd) == Rasters.GRDsource()
@test sourcetrait(:grib) == Rasters.GRIBsource()
@test sourcetrait(:gdal) == Rasters.GDALsource()
@test sourcetrait(:zarr) == Rasters.Zarrsource()

@test sourcetrait(Rasters.NCDsource()) == Rasters.NCDsource()
@test sourcetrait(Rasters.GRDsource()) == Rasters.GRDsource()
@test sourcetrait(Rasters.GRIBsource()) == Rasters.GRIBsource()
@test sourcetrait(Rasters.GDALsource()) == Rasters.GDALsource()
@test sourcetrait(Rasters.Zarrsource()) == Rasters.Zarrsource()

@test sourcetrait(Rasters.NCDsource) == Rasters.NCDsource()
@test sourcetrait(Rasters.GRDsource) == Rasters.GRDsource()
@test sourcetrait(Rasters.GRIBsource) == Rasters.GRIBsource()
@test sourcetrait(Rasters.GDALsource) == Rasters.GDALsource()
@test sourcetrait(Rasters.Zarrsource) == Rasters.Zarrsource()