using Rasters, Zarr
using ZarrDatasets

path = "https://s3.bgc-jena.mpg.de:9000/esdl-esdc-v3.0.2/esdc-16d-2.5deg-46x72x1440-3.0.2.zarr"

@testset "Zarr Raster open" begin
    
r = Raster(path, lazy=true, name="air_temperature_2m")
@test name.(dims(r)) == (:X, :Y, :Ti)
@test length(dims(r, X)) == 144
@test collect(val(dims(r,X))) == collect(-178.75:2.5:178.75)
@test metadata(r)["original_name"] == "t2m"

end
