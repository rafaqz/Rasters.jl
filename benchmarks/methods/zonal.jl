using Rasters, RasterDataSources, ArchGDAL, NaturalEarth

precip = cat((Raster(WorldClim{Climate}, :prec; month = i) for i in 1:12)...; dims = Ti)
all_countries = naturalearth("admin_0_countries", 10)


@be zonal($sum, $(precip[Ti(1)]); of = $countries, progress = $false, threaded = $false) seconds=5
@be zonal($sum, $(precip[Ti(1)]); of = $countries, progress = $false, threaded = $(Rasters.ByMap(false))) seconds=5


@be zonal($sum, $(precip[Ti(1)]); of = $countries, progress = $false, threaded = $true) seconds=5
@be zonal($sum, $(precip[Ti(1)]); of = $countries, progress = $false, threaded = $(Rasters.ByMap(true))) seconds=5