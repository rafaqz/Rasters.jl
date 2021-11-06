using Rasters, Test
using Rasters.LookupArrays, Rasters.Dimensions

@testset "set"  begin
    A = [missing 7; 2 missing]
    ga = Raster(A, (Y(-20:40:20), X(50:10:60)); missingval=missing)

    # Use the Projected lookup, with crs
    lu = Projected(; crs=EPSG(2024))
    ga = set(ga, Y=lu, X=lu)
    @test crs(ga) == EPSG(2024)
    @test mappedcrs(ga) == nothing

    # Set it with usercrs as well
    lu = Projected(; crs=EPSG(2024), mappedcrs=EPSG(4326))
    ga = set(ga, Y=lu, X=lu)
    @test crs(ga) == EPSG(2024)
    @test mappedcrs(ga) == EPSG(4326)

    # Make them intervals
    ga = set(ga, Y=Intervals(Start()), X=Intervals(Start()))
    @test sampling(ga) == (Intervals(Start()), Intervals(Start()))

    # Change the parent array
    ga = set(ga, [1 2; missing 5])
    @test all(parent(ga) .=== [1 2; missing 5])

    # TODO: Add more tests once the lookups are cleaned up
end
