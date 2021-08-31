using GeoData, Test, ArchGDAL

include(joinpath(dirname(pathof(GeoData)), "../test/test_utils.jl"))

A = [missing 7.0f0; 2.0f0 missing]
ga = GeoArray(A, (X, Y); missingval=missing)
ga99 = replace_missing(ga, -9999)
@test all(ga99 .=== [-9999.0f0 7.0f0; 2.0f0 -9999.0f0])
gaNaN = replace_missing(ga, NaN32)
@test all(gaNaN .=== [NaN32 7.0f0; 2.0f0 NaN32])
gaMi = replace_missing(ga, missing)
@test all(gaMi .=== ga)

@testset "boolmask" begin
    @test boolmask(A) == [false true; true false]
    @test boolmask(ga) == [false true; true false]
    @test boolmask(ga99) == [false true; true false]
    @test boolmask(gaNaN) == [false true; true false]
    @test dims(boolmask(ga)) == (X(Base.OneTo(2), mode=NoIndex()), Y(Base.OneTo(2), mode=NoIndex()))
end

@testset "missingmask" begin
    @test all(missingmask(A) .=== [missing true; true missing])
    @test all(missingmask(ga) .=== [missing true; true missing])
    @test all(missingmask(ga99) .=== [missing true; true missing])
    @test all(missingmask(gaNaN) .=== [missing true; true missing])
    @test dims(missingmask(ga)) == (X(Base.OneTo(2), mode=NoIndex()), Y(Base.OneTo(2), mode=NoIndex()))
end

@testset "points" begin
    ga = GeoArray(A, (X(9.0:1.0:10.0), Y(0.1:0.1:0.2)); missingval=missing)
    @test all(points(ga) .=== [(0.1, 9.0) missing; missing (0.2, 10.0)])
    @test all(points(ga; dims=(X, Y)) .=== [(9.0, 0.1) missing; missing (10.0, 0.2)])
    @test all(points(ga; dims=(X, Y), ignore_missing=true) .===
              [(9.0, 0.1) (9.0, 0.2); (10.0, 0.1) (10.0, 0.2)])
end

@testset "trim, crop, extend" begin
    A = [missing missing missing
         missing 2.0     0.5
         missing 1.0     missing]
    ga = GeoArray(A, (X(1.0:1.0:3.0), Y(1.0:1.0:3.0)); missingval=missing)
    # Test with missing on all sides
    ga_r = rot180(ga)
    trimmed = trim(ga)
    trimmed_r = trim(ga_r)
    @test all(trimmed .=== [2.0 0.5; 1.0 missing])
    @test all(trimmed_r .=== [missing 1.0; 0.5 2.0])
    cropped = crop(ga; to=trimmed)
    cropped_r = crop(ga_r; to=trimmed_r)
    @test all(cropped .=== trimmed)
    @test all(cropped_r .=== trimmed_r)
    extended = extend(cropped, ga)[1]
    extended_r = extend(cropped_r; to=ga_r)
    @test all(extended .=== ga)
    @test all(extended_r .=== ga_r)
end


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
    cea = read(geoarray(raster_path))
    GD_output = resample(cea, output_res; crs=output_crs, method=resample_method)

    cea_permuted = permutedims(read(geoarray(raster_path)), (Y, X, Band))
    GD_permuted_output = resample(cea_permuted, output_res; crs=output_crs, method=resample_method)

    # Compare ArchGDAL, resample and permuted resample 
    @test AG_output == GD_output[Band(1)] == permutedims(GD_permuted_output, (X, Y, Band))[Band(1)]
    @test abs(step(dims(GD_output, Y))) ≈ abs(step(dims(GD_output, X))) ≈ 
          abs(step(dims(GD_permuted_output, X))) ≈ output_res

    @testset "snapped size and dim index match" begin
        snaptarget = GD_output
        snapped = resample(cea; to=snaptarget)
        @test size(snapped) == size(snaptarget)
        @test isapprox(index(snaptarget, Y), index(snapped, Y))
        @test isapprox(index(snaptarget, X), index(snapped, X))
    end
end
