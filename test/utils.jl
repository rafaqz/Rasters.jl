using Rasters, Test
using Rasters: cleankeys

@testset "cleankeys" begin
    @test cleankeys(["a", "b", "c"]) == (:a, :b, :c)
    @test cleankeys(("a", "b", "c")) == (:a, :b, :c)
end
