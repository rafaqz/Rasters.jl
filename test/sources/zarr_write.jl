using Test
using Rasters
using DimensionalData
using ZarrDatasets

@testset "zarr write with force" begin
    ras = Raster(rand(X(0:0.1:1), Y(0:0.2:1)))

    mktempdir() do tmpdir
        fn = joinpath(tmpdir, "test.zarr")

        # Initial write
        write(fn, ras)
        @test isdir(fn)
        @test all(Raster(fn) .=== ras)

        # Without force should error
        @test_throws ArgumentError write(fn, ras)

        # force=true should overwrite (issue #1060)
        write(fn, ras; force=true)
        @test all(Raster(fn) .=== ras)

        # Overwrite with different data
        ras2 = rebuild(ras; data=ones(size(ras)))
        write(fn, ras2; force=true)
        @test all(Raster(fn) .=== ras2)
    end
end
