using Rasters, Test, ArchGDAL, Dates, Statistics
using Rasters.LookupArrays, Rasters.Dimensions

include(joinpath(dirname(pathof(Rasters)), "../test/test_utils.jl"))

A = [missing 7.0f0; 2.0f0 missing]
ga = Raster(A, (X, Y); missingval=missing)
ga99 = replace_missing(ga, -9999)
@test all(ga99 .=== [-9999.0f0 7.0f0; 2.0f0 -9999.0f0])
gaNaN = replace_missing(ga, NaN32)
@test all(gaNaN .=== [NaN32 7.0f0; 2.0f0 NaN32])
gaMi = replace_missing(ga, missing)
@test all(gaMi .=== ga)

@testset "boolmask" begin
    @test boolmask(Array, ga) == [false true; true false]
    @test boolmask(ga) == [false true; true false]
    @test boolmask(ga99) == [false true; true false]
    @test boolmask(gaNaN) == [false true; true false]
    @test dims(boolmask(ga)) == (X(NoLookup(Base.OneTo(2))), Y(NoLookup(Base.OneTo(2))))
end

@testset "missingmask" begin
    @test all(missingmask(A) .=== [missing true; true missing])
    @test all(missingmask(ga) .=== [missing true; true missing])
    @test all(missingmask(ga99) .=== [missing true; true missing])
    @test all(missingmask(gaNaN) .=== [missing true; true missing])
    @test dims(missingmask(ga)) == (X(NoLookup(Base.OneTo(2))), Y(NoLookup(Base.OneTo(2))))
end

@testset "mask" begin
    A1 = [missing 1; 2 3]
    ga1 = Raster(A1, (X, Y); missingval=missing)
    @test all(mask(ga1; to=ga) .=== [missing 1; 2 missing])
    ga2 = replace_missing(ga1 .* 1.0; missingval=NaN)
    @test all(mask(ga2; to=ga) .=== [NaN 1.0; 2.0 NaN])
    ga3 = replace_missing(ga1; missingval=-9999)
    mask!(ga3; to=ga)
    @test all(ga3 .=== [-9999 1; 2 -9999])
end

@testset "classify" begin
    A1 = [missing 1; 2 3]
    ga1 = Raster(A1, (X, Y); missingval=missing)
    @test all(classify(ga1, 1=>99, 2=>88, 3=>77) .=== [missing 99; 88 77])
    @test all(classify(ga1, 1=>99, 2=>88, 3=>77; others=0) .=== [missing 99; 88 77])
    @test all(classify(ga1, 1=>99, 2=>88; others=0) .=== [missing 99; 88 0])

    A2 = [1.0 2.5; 3.0 4.0]
    ga2 = Raster(A2, (X, Y); missingval=missing)
    @test classify(ga2, (2, 3)=>:x, >(3)=>:y) == [1.0 :x; 3.0 :y]
    @test classify(ga2, (>=(1), <(2))=>:x, >=(3)=>:y) == [:x 2.5; :y :y]
    classify!(ga2, (1, 2.5)=>0.0, >=(3)=>-1.0; lower=(>), upper=(<=))
    @test ga2 == [1.0 0.0; -1.0 -1.0]
    classify!(ga2, [1 2.5 0.0; 2.5 4.0 -1.0]; lower=(>), upper=(<=))
    @test ga2 == [1.0 0.0; -1.0 -1.0]
end

@testset "points" begin
    ga = Raster(A, (X(9.0:1.0:10.0), Y(0.1:0.1:0.2)); missingval=missing)
    @test all(collect(points(ga; order=(Y, X))) .=== [missing (0.2, 9.0); (0.1, 10.0) missing])
    @test all(collect(points(ga; order=(X, Y))) .=== [missing (9.0, 0.2); (10.0, 0.1) missing])
    @test all(points(ga; order=(X, Y), ignore_missing=true) .===
              [(9.0, 0.1) (9.0, 0.2); (10.0, 0.1) (10.0, 0.2)])
end

@testset "extract" begin
    A = [1 2; 3 4]
    ga = Raster(A, (X(9.0:1.0:10.0), Y(0.1:0.1:0.2)); name=:test, missingval=missing)
    @testset "points" begin
        @test all(extract(ga, [(9.0, 0.1), (10.0, 0.2), (10.0, 0.3)]) .=== 
                  [(X=9.0, Y=0.1, test=1), (X=10.0, Y=0.2, test=4), (X=10.0, Y=0.3, test=missing)])
        @test all(extract(ga, [(0.1, 9.0), (0.2, 10.0), (0.3, 10.0)]; order=(Y, X)) .=== 
                  [(Y=0.1, X=9.0, test=1), (Y=0.2, X=10.0, test=4), (Y=0.3, X=10.0, test=missing)])
    end
    @testset "Tables.jl compatible" begin
        @test all(extract(ga, [(X=9.0, Y=0.1), (X=10.0, Y=0.2), (X=10.0, Y=0.3)]) .=== 
              [(X=9.0, Y=0.1, test=1), (X=10.0, Y=0.2, test=4), (X=10.0, Y=0.3, test=missing)])
        @test all(extract(ga, [(X=9.0, Y=0.1), (X=10.0, Y=0.2), (X=10.0, Y=0.3)]; order=(Y=>:Y, X=>:X)) .=== 
              [(Y=0.1, X=9.0, test=1), (Y=0.2, X=10.0, test=4), (Y=0.3, X=10.0, test=missing)])
    end
end

@testset "trim, crop, extend" begin
    A = [missing missing missing
         missing 2.0     0.5
         missing 1.0     missing]
    ga = Raster(A, (X(1.0:1.0:3.0), Y(1.0:1.0:3.0)); missingval=missing)
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

@testset "mosaic" begin
    reg1 = Raster([0.1 0.2; 0.3 0.4], (X(2.0:1.0:3.0), Y(5.0:1.0:6.0)))
    reg2 = Raster([1.1 1.2; 1.3 1.4], (X(3.0:1.0:4.0), Y(6.0:1.0:7.0)))
    irreg1 = Raster([0.1 0.2; 0.3 0.4], (X([2.0, 3.0]), Y([5.0, 6.0])))
    irreg2 = Raster([1.1 1.2; 1.3 1.4], (X([3.0, 4.0]), Y([6.0, 7.0])))

    span_x1 = Explicit(vcat((1.5:1.0:2.5)', (2.5:1.0:3.5)'))
    span_x2 = Explicit(vcat((2.5:1.0:3.5)', (3.5:1.0:4.5)'))
    exp1 = Raster([0.1 0.2; 0.3 0.4], (X(Sampled([2.0, 3.0]; span=span_x1)), Y([5.0, 6.0])))
    exp2 = Raster([1.1 1.2; 1.3 1.4], (X(Sampled([3.0, 4.0]; span=span_x2)), Y([6.0, 7.0])))
    @test val(span(mosaic(first, exp1, exp2), X)) == [1.5 2.5 3.5; 2.5 3.5 4.5]
    @test all(mosaic(first, reg1, reg2) .=== 
              mosaic(first, irreg1, irreg2) .===
              mosaic(first, irreg1, irreg2) .=== 
              [0.1 0.2 missing; 
               0.3 0.4 1.2; 
               missing 1.3 1.4])
    @test all(mosaic(last, reg1, reg2) .===
              mosaic(last, irreg1, irreg2) .===
              mosaic(last, exp1, exp2) .=== [0.1 0.2 missing; 
                                             0.3 1.1 1.2; 
                                             missing 1.3 1.4])

    # 3 dimensions
    A1 = Raster(ones(2, 2, 2), (X(2.0:-1.0:1.0), Y(5.0:1.0:6.0), Ti(DateTime(2001):Year(1):DateTime(2002))))
    A2 = Raster(zeros(2, 2, 2), (X(3.0:-1.0:2.0), Y(4.0:1.0:5.0), Ti(DateTime(2002):Year(1):DateTime(2003))))
    @test all(mosaic(mean, A1, A2) |> parent .=== cat([missing missing missing
                                                 missing 1.0     1.0
                                                 missing 1.0     1.0    ],
                                                [0.0     0.0 missing
                                                 0.0     0.5     1.0   
                                                 missing 1.0     1.0    ],
                                                [0.0     0.0     missing
                                                 0.0     0.0     missing    
                                                 missing missing missing], dims=3))
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
