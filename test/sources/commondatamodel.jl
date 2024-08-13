using Rasters, NCDatasets, Test
import Rasters: ForwardOrdered, ReverseOrdered, Regular
@testset "step" begin
    # test if regular indices are correctly rounded
    f32_indices = range(0.075f0, 10.075f0; step = 0.05f0) |> collect
    @test Rasters._cdmspan(f32_indices, ForwardOrdered())[1] === Regular(0.05)

    f32_indices_rev = range(10.075f0, 0.075f0; step = -0.05f0) |> collect
    @test Rasters._cdmspan(f32_indices_rev, ReverseOrdered())[1] === Regular(-0.05)

    # test if regular indices are not rounded when they should not
    indices_one_third = range(0, 10; length = 31) |> collect
    @test Rasters._cdmspan(indices_one_third, ForwardOrdered())[1] === Regular(1/3)

    # test when reading a file
    ras = Raster(rand(X(f32_indices), Y(indices_one_third)))
    tempfile = tempname() * ".nc"
    write(tempfile, ras)
    ras_read = Raster(tempfile)
    steps = step.(dims(ras_read))
    @test steps[1] == 0.05
    @test steps[2] == 1/3

end