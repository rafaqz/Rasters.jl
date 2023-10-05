using Rasters, Rasters.LookupArrays, ArchGDAL
using Test

include(joinpath(dirname(pathof(Rasters)), "../test/test_utils.jl"))

@testset "cellsize" begin
    dimz = X(Projected(90.0:0.1:99.9; sampling=Intervals(Start()), order = ForwardOrdered(), span = Regular(0.1), crs=EPSG(4326))),
       Y(Projected(0.0:0.1:89.9; sampling=Intervals(Start()), order = ForwardOrdered(), span = Regular(0.1), crs=EPSG(4326)))

    dimz_25832 = X(Projected(0.0:100:10000.0; sampling=Intervals(Start()), order = ForwardOrdered(), span = Regular(100), crs=EPSG(25832))),
       Y(Projected(0.0:100:10000.0; sampling=Intervals(Start()), order = ForwardOrdered(), span = Regular(100), crs=EPSG(25832)))

    cs = cellsize(dimz)
    cs2 = cellsize(dimz_25832)

    # Check the output is a raster 
    @test cs isa Raster
    @test cs2 isa Raster
    # Test that the total area matches the expected area (1/72th of the Earth surface)
    @test sum(cs) ≈ 510.1e6/72 rtol = 0.01
    # Test all areas are about 0.01 km2
    @test maximum(cs2) ≈ 0.01 rtol = 0.01
    @test minimum(cs2) ≈ 0.01 rtol = 0.01

end
