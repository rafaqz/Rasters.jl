using GeoData, Test

A = [missing 7.0f0; 2.0f0 missing]
ga = GeoArray(A, (X, Y); missingval=missing)
ga99 = replace_missing(ga, -9999)
@test all(ga99 .=== [-9999.0f0 7.0f0; 2.0f0 -9999.0f0])
gaNaN = replace_missing(ga, NaN32)
@test all(gaNaN .=== [NaN32 7.0f0; 2.0f0 NaN32])

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
