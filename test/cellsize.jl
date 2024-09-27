using Rasters, DimensionalData, Rasters.Lookups, ArchGDAL
using Test
using DimensionalData: @dim, YDim
include(joinpath(dirname(pathof(Rasters)), "../test/test_utils.jl"))

@testset "cellsize" begin
    x = X(Projected(90.0:0.1:99.9; sampling=Intervals(Start()), order = ForwardOrdered(), span = Regular(0.1), crs=EPSG(4326)))
    y = Y(Projected(0.0:0.1:89.9; sampling=Intervals(Start()), order = ForwardOrdered(), span = Regular(0.1), crs=EPSG(4326)))
    y_rev = Y(Projected(89.9:-0.1:0; sampling=Intervals(Start()), order = ReverseOrdered(), span = Regular(-0.1), crs=EPSG(4326)))
    dimz = (x,y)

    dimz_25832 = X(Projected(0.0:100:10000.0; sampling=Intervals(Start()), order = ForwardOrdered(), span = Regular(100), crs=EPSG(25832))),
       Y(Projected(0.0:100:10000.0; sampling=Intervals(Start()), order = ForwardOrdered(), span = Regular(100), crs=EPSG(25832)))
    ras = ones(dimz)

    cs = cellsize(dimz)
    cs2 = cellsize(dimz_25832)
    cs3 = cellsize((x, y_rev))

    # Check the output is a raster 
    @test cs isa Raster
    @test cs2 isa Raster
    # Test that the total area matches the expected area (1/72th of the Earth surface)
    @test sum(cs) ≈ sum(cs3)
    @test sum(cs) ≈ 510.1e6/72 rtol = 0.01
    # Test all areas are about 0.01 km2
    @test maximum(cs2) ≈ 0.01 rtol = 0.01
    @test minimum(cs2) ≈ 0.01 rtol = 0.01
    # test passing in a raster or dims gives the same result
    cs_ras = cellsize(ras)
    @test cs == cs_ras

    # test a YDim other than Y is handled correctly
    @dim Lat YDim "latitude"
    lat = Lat(Projected(0.0:0.1:89.9; sampling=Intervals(Start()), order = ForwardOrdered(), span = Regular(0.1), crs=EPSG(4326)))
    cs_lat = cellsize((x, lat))
    @test parent(cs_lat) == parent(cs)
    @test dims(cs_lat) == (x, lat)

    # test point sampling throws an error
    pointsy = set(y, Points())
    @test_throws ArgumentError cellsize((x, pointsy))

    # test missing crs throws an error
    nocrsdimz = setcrs(dimz, nothing)
    @test_throws ArgumentError cellsize(nocrsdimz)

end