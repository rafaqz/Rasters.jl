using GeoData, Test

@testset "set"  begin
    A = [missing 7; 2 missing]
    ga = GeoArray(A, (Lat(-20:40:20), Lon(50:10:60)); missingval=missing)

    # Use the Projected mode, with crs
    mode = Projected(crs=EPSG(2024))
    ga = set(ga, Lat=mode, Lon=mode)
    @test crs(ga) == EPSG(2024)
    @test mappedcrs(ga) == nothing

    # Set it with usercrs as well
    mode = Projected(crs=EPSG(2024), mappedcrs=EPSG(4326))
    ga = set(ga, Lat=mode, Lon=mode)
    @test crs(ga) == EPSG(2024)
    @test mappedcrs(ga) == EPSG(4326)

    # Make them intervals
    ga = set(ga, Lat=Intervals(Start()), Lon=Intervals(Start()))
    @test sampling(ga) == (Intervals(Start()), Intervals(Start()))

    # Change the parent array
    ga = set(ga, [1 2; missing 5])
    @test all(parent(ga) .=== [1 2; missing 5])

    # TODO: Add more tests once the modes are cleaned up
end
