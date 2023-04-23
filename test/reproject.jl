using Rasters, Test, ArchGDAL
using Rasters.LookupArrays, Rasters.Dimensions
using Rasters: reproject, convertlookup

@testset "reproject" begin
    cea = ProjString("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
    wktcea = convert(WellKnownText, cea)
    projcea = convert(ProjString, cea)
    wkt4326 = convert(WellKnownText, EPSG(4326))
    proj4326 = convert(ProjString, EPSG(4326))

    @test reproject(proj4326, cea, Y(), -30.0) ≈ -3.658789324855012e6
    @test reproject(wktcea, EPSG(4326), Y(), -3.658789324855012e6) ≈ -30.0
    @test reproject(projcea, wkt4326, Y(), -3.658789324855012e6) ≈ -30.0
    @test reproject(cea, proj4326, Y(), [-3.658789324855012e6]) ≈ [-30.0]

    @test reproject(proj4326, cea, X(), 180.0) ≈ 1.7367530445161372e7
    @test reproject(cea, EPSG(4326), X(), 1.7367530445161372e7) ≈ 180.0
    @test reproject(cea, wkt4326, X(), 1.7367530445161372e7) ≈ 180.0
    @test reproject(cea, proj4326, X(), [1.7367530445161372e7]) ≈ [180.0]

    x = X(Projected(0.0:1.0:360.0; crs=EPSG(4326), dim=X(), order=ForwardOrdered(), span=Regular(1.0), sampling=Intervals(Start())))
    y = Y(Projected(-80.0:1.0:80.0; crs=EPSG(4326), dim=Y(), order=ReverseOrdered(), span=Regular(1.0), sampling=Intervals(Start())))
    x1, y1 = reproject(projcea, (x, y))
    @test span(x1) isa Irregular
    @test span(y1) isa Irregular
    x2, y2 = reproject(EPSG(4326), (x, y))
    @test all(x .== x2)
    @test all(y .== y2)
    @test span(x2) == Regular(1.0)
    @test span(y2) == Regular(1.0)
end

@testset "convertlookup" begin
    projcea = ProjString("+proj=cea")
    proj4326 = convert(ProjString, EPSG(4326))

    lonstart, lonend = 0.5, 179.5
    cealonstart, cealonend = reproject(proj4326, projcea, X(), [lonstart, lonend])
    cealonrange = cealonstart:1.0:cealonend
    lon = X(Projected(cealonrange, ForwardOrdered(), Regular(step(cealonrange)), 
                      Intervals(Center()), NoMetadata(), projcea, proj4326, X()))
    convertedlon = convertlookup(Mapped, lon)
    @test all(isapprox.(bounds(convertedlon), (0.0, 180.0); atol=1e-10))
    @test val(convertedlon) ≈ 0.5:179.5

    projectedlon = convertlookup(Projected, convertedlon)
    @test val(projectedlon) ≈ val(lon)
    @test all(isapprox.(bounds(projectedlon), bounds(lon); atol=1e-10))

    # Y in not linear
    latbounds = -90.0, 90.0
    ceabounds = reproject(proj4326, projcea, Y(), latbounds)
    # Vals are the bounding range without the end bound
    cealatrange = ceabounds[1] : (ceabounds[2] - ceabounds[1] / 180) : ceabounds[2]
    lat = Y(Projected(cealatrange, ForwardOrdered(), Regular(step(cealatrange)), 
                      Intervals(Start()), NoMetadata(), projcea, proj4326, Y()))
    convertedlat = convertlookup(Mapped, lat)
    @test first(convertedlat) ≈ -90.0 atol=1e-5
    @test all(isapprox.(bounds(convertedlat), (-90.0, 90.0); atol=1e-5))
    @test val(lat) ≈ reproject(proj4326, projcea, Y(), val(convertedlat))

    projectedlat = convertlookup(Projected, convertedlat)
    @test val(projectedlat) ≈ val(lat)
    @test all(bounds(projectedlat) .≈ bounds(lat))

    A = DimArray(zeros(length(lon), length(lat)), (lon, lat))
    Aconv = convertlookup(Mapped, A)

    @test index(Aconv) == (index(convertedlon), index(convertedlat))
    @test val.(span(Aconv)) == val.(span.((convertedlon, convertedlat)))
end
