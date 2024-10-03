using Rasters, DimensionalData, Rasters.Lookups, ArchGDAL
using Test
using DimensionalData: @dim, YDim
include(joinpath(dirname(pathof(Rasters)), "../test/test_utils.jl"))

@testset "cellarea" begin
    x = X(Projected(90.0:0.1:99.9; sampling=Intervals(Start()), order = ForwardOrdered(), span = Regular(0.1), crs=EPSG(4326), dim = X()))
    y = Y(Projected(0.0:0.1:89.9; sampling=Intervals(Start()), order = ForwardOrdered(), span = Regular(0.1), crs=EPSG(4326), dim = Y()))
    y_rev = Y(Projected(89.9:-0.1:0; sampling=Intervals(Start()), order = ReverseOrdered(), span = Regular(-0.1), crs=EPSG(4326), dim = Y()))
    dimz = (x,y)
    ras = ones(dimz)

    dimz_25832 = X(Projected(0.0:100:10000.0; sampling=Intervals(Start()), order = ForwardOrdered(), span = Regular(100), crs=EPSG(25832))),
       Y(Projected(0.0:100:10000.0; sampling=Intervals(Start()), order = ForwardOrdered(), span = Regular(100), crs=EPSG(25832)))
    # to test that this is still stable on as 5x5m grid
    smallspan_25832 = X(Projected(0.0:5:100.0; sampling=Intervals(Start()), order = ForwardOrdered(), span = Regular(5), crs=EPSG(25832))),
       Y(Projected(0.0:5:100.0; sampling=Intervals(Start()), order = ForwardOrdered(), span = Regular(5), crs=EPSG(25832)))

    cs = cellarea(dimz)
    cs2 = cellarea(dimz_25832)
    cs3 = cellarea((x, y_rev))
    cs4 = cellarea(smallspan_25832)

    # check cellsize still works
    csold = @test_warn "deprecated" cellsize(dimz)
    @test all(csold .≈ (cs ./ 1e6))

    # Check the output is a raster 
    @test cs isa Raster
    @test cs2 isa Raster
    # Test that the total area matches the expected area (1/72th of the Earth surface)
    @test sum(cs) ≈ sum(cs3)
    @test sum(cs) ≈ 510.1e12/72 rtol = 0.01
    # Test all areas are about 10,000 m2 for the 100x100m grid
    @test maximum(cs2) ≈ 1e4 rtol = 0.01
    @test minimum(cs2) ≈ 1e4 rtol = 0.01
    # Test all areas are about 25 m2 on the 5x5m grid
    @test maximum(cs4) ≈ 25 rtol = 0.01
    @test minimum(cs4) ≈ 25 rtol = 0.01
    # test passing in a raster or dims gives the same result
    cs_ras = cellarea(ras)
    @test cs == cs_ras

    ## test with area_in_crs true
    cs_planar = cellarea(dimz; area_in_crs = true)
    cs_planar2 = cellarea(dimz_25832; area_in_crs = true)
    cs_planar3 = cellarea((x, y_rev); area_in_crs = true)
    @test all(≈(0.01), cs_planar)
    @test all(≈(10_000), cs_planar2)
    @test all(≈(0.01), cs_planar3)

    # test a YDim other than Y is handled correctly
    @dim Lat YDim "latitude"
    lat = Lat(Projected(0.0:0.1:89.9; sampling=Intervals(Start()), order = ForwardOrdered(), span = Regular(0.1), crs=EPSG(4326)))
    cs_lat = cellarea((x, lat))
    @test parent(cs_lat) == parent(cs)
    @test dims(cs_lat) == (x, lat)

    # test using mapped dims to calculate cellsize
    x_3857 = X(Projected(0:1e5:1e7; sampling=Intervals(Start()), order = ForwardOrdered(), span = Regular(1e5), crs=EPSG(3857), dim = X()))
    y_3857 = Y(Projected(0:1e5:1e7; sampling=Intervals(Start()), order = ForwardOrdered(), span = Regular(1e5), crs=EPSG(3857), dim = Y()))
    dimz_3857 = (x_3857, y_3857)
    mappeddimz = setmappedcrs(dimz_3857, EPSG(4326))
    @test all(isapprox.(cellarea((dimz_3857)), cellarea((mappeddimz)); rtol = 0.001))
    @test parent(cellarea((dimz_3857))) != parent(cellarea((mappeddimz))) # to make sure it actually used the mapped dims
    # test point sampling throws an error
    pointsy = set(y, Points())
    @test_throws ArgumentError cellarea((x, pointsy))

    # test missing crs throws an error
    nocrsdimz = setcrs(dimz, nothing)
    @test_throws ArgumentError cellarea(nocrsdimz)
end