using GeoData, Test, ArchGDAL
using GeoData: reproject, convertmode

@testset "reproject" begin
    cea = ProjString("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
    wktcea = convert(WellKnownText, cea)
    projcea = convert(ProjString, cea)
    wkt4326 = convert(WellKnownText, EPSG(4326))
    proj4326 = convert(ProjString, EPSG(4326))

    @test reproject(proj4326, cea, Lat(), -30.0) ≈ -3.658789324855012e6
    @test reproject(wktcea, EPSG(4326), Lat(), -3.658789324855012e6) ≈ -30.0
    @test reproject(projcea, wkt4326, Lat(), -3.658789324855012e6) .≈ -30.0
    @test reproject(cea, proj4326, Lat(), [-3.658789324855012e6]) ≈ [-30.0]

    @test reproject(proj4326, cea, Lon(), 180.0) ≈ 1.7367530445161372e7
    @test reproject(cea, EPSG(4326), Lon(), 1.7367530445161372e7) ≈ 180.0
    @test reproject(cea, wkt4326, Lon(), 1.7367530445161372e7) .≈ 180.0
    @test reproject(cea, proj4326, Lon(), [1.7367530445161372e7]) ≈ [180.0]
end

@testset "convertmode" begin
    projcea = ProjString("+proj=cea")
    proj4326 = convert(ProjString, EPSG(4326))

    lonstart, lonend = 0.5, 179.5
    cealonstart, cealonend = reproject(proj4326, projcea, Lon(), [lonstart, lonend])
    cealonrange = LinRange(cealonstart, cealonend, 180)
    lon = Lon(cealonrange; mode=Projected(Ordered(), Regular(step(cealonrange)), 
              Intervals(Center()), projcea, proj4326))
    convertedlon = convertmode(Mapped, lon)
    @test all(isapprox.(bounds(convertedlon), (0.0, 180.0); atol=1e-10))
    @test val(convertedlon) ≈ 0.5:179.5

    projectedlon = convertmode(Projected, convertedlon)
    @test val(projectedlon) ≈ val(lon)
    @test all(isapprox.(bounds(projectedlon), bounds(lon); atol=1e-10))

    # Lat in not linear
    latbounds = -90.0, 90.0
    ceabounds = reproject(proj4326, projcea, Lat(), latbounds)
    # Vals are the bounding range without the end bound
    cealatrange = LinRange(ceabounds..., 181)[1:180]
    lat = Lat(cealatrange; mode=Projected(Ordered(), Regular(step(cealatrange)), 
              Intervals(Start()), projcea, proj4326))
    convertedlat = convertmode(Mapped, lat)
    @test first(convertedlat) ≈ -90.0 atol=1e-5
    @test all(isapprox.(bounds(convertedlat), (-90.0, 90.0); atol=1e-5))
    @test val(lat) ≈ reproject(proj4326, projcea, Lat(), val(convertedlat))

    projectedlat = convertmode(Projected, convertedlat)
    @test val(projectedlat) ≈ val(lat)
    @test all(bounds(projectedlat) .≈ bounds(lat))

    A = DimArray(zeros(length(lon), length(lat)), (lon, lat))
    Aconv = convertmode(Mapped, A)

    @test dims(Aconv) == (convertedlon, convertedlat)
end
