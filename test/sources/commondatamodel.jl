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
    write(tempfile, ras; force=true)
    ras_read = Raster(tempfile)
    steps = step.(dims(ras_read))
    @test steps[1] == 0.05
    @test steps[2] == 1/3

end

@testset "faulty / no grid mapping" begin
    # Construct a NetCDF dataset
    filepath = joinpath(tempdir(), "test.nc")
    ds = NCDataset(filepath, "c")
    data = [Float32(i + j) for i = 1:100, j = 1:110]
    v = defVar(ds, "temperature", data, ("lon", "lat"))
    # Give `v` a "grid_mapping" attribute
    # that points to a non-existent variable
    v.attrib["grid_mapping"] = "non_existent_variable"
    close(ds)
    
    # try loading raster
    @test_nowarn Raster(filepath)
    # make sure we don't magically materialize a CRS
    @test isnothing(Rasters.crs(Raster(filepath)))

    # Once Rasters has better CF-CRS support,
    # we should be able to load a splatted global CRS,
    # or a CRS as a global attribute a la Zarr.
end
