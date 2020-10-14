using GeoData, Test, ArchGDAL
using GeoFormatTypes
using GeoData: resample

@testset "resample" begin
    download("https://download.osgeo.org/geotiff/samples/gdal_eg/cea.tif", "data/cea.tif")

    raster_path = "data/cea.tif"
    output_res = 0.0027
    output_crs = EPSG(4326)
    resample_method = "near"

    ## Resample cea.tif manually with ArchGDAL
    wkt = convert(String, convert(WellKnownText, output_crs))
    AG_output = ArchGDAL.read(raster_path) do dataset
        ArchGDAL.gdalwarp([dataset], ["-t_srs", "$(wkt)",
                                "-tr", "$(output_res)", "$(output_res)",
                                "-r", "$(resample_method)"]) do warped
            ArchGDAL.read(ArchGDAL.getband(warped, 1))
        end
    end

    ## Resample cea.tif using resample
    cea = GeoArray(GDALarray(raster_path))
    GD_output = resample(cea, output_res, crs = output_crs, method = resample_method)

    ## Compare the two
    @test AG_output == GD_output.data[:, :, 1]
    @test abs(step(dims(GD_output)[1])) ≈ abs(step(dims(GD_output)[2])) ≈ output_res
end