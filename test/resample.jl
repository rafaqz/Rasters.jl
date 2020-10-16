using GeoData, Test, ArchGDAL, GeoFormatTypes
using GeoData: resample

include(joinpath(dirname(pathof(GeoData)), "../test/test_utils.jl"))

@testset "resample" begin
    raster_path = maybedownload("https://download.osgeo.org/geotiff/samples/gdal_eg/cea.tif")

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
    projectedcrs(cea)
    GD_output = resample(cea, output_res, crs=output_crs, method=resample_method)

    ## Compare the two
    @test AG_output == GD_output.data[:, :, 1]
    @test abs(step(dims(GD_output, Lat))) ≈ abs(step(dims(GD_output, Lon))) ≈ output_res

    @testset "snapped size and dim index match" begin
        snaptarget = GD_output
        snapped = resample(cea, snaptarget)
        @test size(snapped) == size(snaptarget)
        @test isapprox(index(snaptarget, Lat), index(snapped, Lat))
        @test isapprox(index(snaptarget, Lon), index(snapped, Lon))
    end
end
