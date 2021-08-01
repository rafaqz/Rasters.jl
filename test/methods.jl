using GeoData, Test

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
    extended = extend(cropped; to=ga)
    extended_r = extend(cropped_r; to=ga_r)
    @test all(extended .=== ga)
    @test all(extended_r .=== ga_r)
end
