using Rasters, Test, Dates, DimensionalData
using Rasters: FileArray, ASCIIfile

const DD = DimensionalData

ascpath = "britanny.asc"

@testset "ASCII array" begin
    @time ascarray = Raster(ascpath)
    @time lazyarray = Raster(ascpath; lazy=true)
    @time eagerarray = Raster(ascpath; lazy=false)
    @test_throws ArgumentError Raster("notafile.asc")

    @testset "lazyness" begin
        @test parent(ascarray) isa Array
        @test parent(lazyarray) isa FileArray
        @test parent(eagerarray) isa Array
    end

    @testset "open" begin
        @test all(open(A -> A[Y=1], ascarray) .=== ascarray[:, 1, :])
        tempfile = tempname()
        cp(ascpath, tempfile * ".asc")
        ascwritearray = Raster(tempfile * ".asc"; lazy=true)
        open(ascwritearray; write=true) do A
            A .*= 2
        end
        @test Raster(tempfile * ".asc") == ascarray .* 2
    end

    @testset "read" begin
        A = read(ascarray)
        @test A isa Raster
        @test parent(A) isa Array
        A2 = zero(A)
        @time read!(ascarray, A2);
        A3 = zero(A)
        @time read!(ascpath, A3);
        @test A == A2 == A3
    end
end