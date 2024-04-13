using Rasters, Test
using Rasters: cleankeys, _chunks_to_tuple
using Rasters.DiskArrays

@testset "cleankeys" begin
    @test cleankeys(["a", "b", "c"]) == (:a, :b, :c)
    @test cleankeys(("a", "b", "c")) == (:a, :b, :c)
end

@testset "_chunks_to_tuple" begin
    A = zeros(100, 100)
    template = DiskArrays.RechunkedDiskArray(A, DiskArrays.GridChunks(A, (32, 16, 1)))
    _chunks_to_tuple(template, (X(), Y(), Ti()), true) == (32, 16, 1)
    @test _chunks_to_tuple(template, (X(), Y()), true) == (32, 16)
    @test _chunks_to_tuple(template, (Y(), X()), true) == (32, 16)
    @test _chunks_to_tuple(template, (Y(), X()), false) == nothing

    @test_throws ArgumentError _chunks_to_tuple(template, (Y(), X()), (16, 32, 1))
    @test_throws ArgumentError _chunks_to_tuple(template, (Y(), X()), (16,))
    @test _chunks_to_tuple(template, (Y(), X()), (16, 32)) == (16, 32)

    @test _chunks_to_tuple(template, (X(), Y()), (X=16, Y=32)) == (16, 32)
    @test _chunks_to_tuple(template, (X(), Y()), (Y=16, X=32)) == (32, 16)

    @test _chunks_to_tuple(template, (X(), Y()), (Y(8), X(12))) == (12, 8)
    @test _chunks_to_tuple(template, (X(), Y()), (Y(8), X(12), Ti(5))) == (12, 8)
    @test _chunks_to_tuple(template, (X(), Y(), Ti()), (Y(8), X(12))) == (12, 8, 1)
    @test _chunks_to_tuple(template, (Ti(), X(), Y()), (Y(8), X(12))) == (1, 12, 8)
end
