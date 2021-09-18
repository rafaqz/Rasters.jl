using GeoData, RasterDataSources, Test, Dates, NCDatasets, ArchGDAL

# Too big to test on CI
# if !haskey(ENV, "CI")
#     @testset "load WorldClim Weather" begin
#         # Weather time-series
#         dates = (Date(2001), Date(2002))
#         ser = series(WorldClim{Weather}, (:prec,); date=dates)
#         ser[Date(2001, 1)][:prec]
#         A = geoarray(WorldClim{Weather}, :prec; date=DateTime(2001, 05), mappedcrs=EPSG(4326))
#     end
# end

@testset "load WorldClim Climate" begin
    # Weather time-series
    ser = GeoSeries(WorldClim{Climate}, :prec; res="10m", month=Jan:March, mappedcrs=EPSG(4326))
    # Select Australia, using regular lat/lon selectors
    A = ser[month=Jan]
    @test A isa GeoArray
    A[Y(Between(-10, -45)), X(Between(110, 160))]
    st = GeoStack(WorldClim{Climate}, (:prec, :tmax); month=1)
    st[:prec]
    @test st isa GeoStack
end

@testset "load WorldClim BioClim" begin
    A = GeoArray(WorldClim{BioClim}, :Bio_1; mappedcrs=EPSG(4326))
    A[Y(Between(-10, -45)), X(Between(110, 160))]
    @test A isa GeoArray
    st = stack(WorldClim{BioClim}, (1, 2))
    st[:bio1]
    @test st isa GeoStack
    @test A isa GeoArray
end

@testset "load CHELSA BioClim" begin
    A = GeoArray(CHELSA{BioClim}, 1; mappedcrs=EPSG(4326))
    @test GeoData.name(A) == :bio1
    st = GeoStack(CHELSA{BioClim}, (:bio1, :BIO2))
    @test keys(st) == (:bio1, :bio2)
    @test A isa GeoArray
    @test st isa GeoStack
    @test st[:bio2] isa GeoArray
end

@testset "load EarthEnv HabitatHeterogeneity" begin
    A = GeoArray(EarthEnv{HabitatHeterogeneity}, :cv; mappedcrs=EPSG(4326))
    A[Y(Between(-10, -45)), X(Between(110, 160))] 
    st = GeoStack(EarthEnv{HabitatHeterogeneity}, (:cv, :evenness))
    @test A isa GeoArray
    @test st isa GeoStack
    @test st[:evenness] isa GeoArray
end

@testset "load EarthEnv LandCover" begin
    A = GeoArray(EarthEnv{LandCover}, 2; mappedcrs=EPSG(4326))
    @test GeoData.name(A) == :evergreen_broadleaf_trees
    A[Y(Between(-10, -45)), X(Between(110, 160))]
    @test A isa GeoArray
    st = GeoStack(EarthEnv{LandCover}, (:evergreen_broadleaf_trees, :deciduous_broadleaf_trees); mappedcrs=EPSG(4326))
    keys(st) = (:evergreen_broadleaf_trees, :deciduous_broadleaf_trees)
    @test st isa GeoStack
end

@testset "load ALWB" begin
    A = GeoArray(ALWB{Deciles,Day}, :rain_day; date=DateTime(2019, 10, 19))
    @test crs(A) == EPSG(4326)
    A = GeoArray(ALWB{Values,Day}, :ss_pct; date=DateTime(2019, 10, 19))
    @test crs(A) == EPSG(4326)
    st = GeoStack(ALWB{Values,Day}, (:s0_pct, :ss_pct); date=DateTime(2019, 10, 19))
    @test crs(st) == EPSG(4326)
    @test crs(st[:s0_pct]) == EPSG(4326)
    dates = DateTime(2019, 10, 19), DateTime(2021, 11, 20)
    s = GeoSeries(ALWB{Values,Day}, (:s0_pct, :ss_pct); date=dates)
    s[1]
    @test A isa GeoArray
    @test st isa GeoStack
    @test s isa GeoSeries
end

# Obscure .Z format may not work on windows
if Sys.islinux()
    @testset "load AWAP" begin
        A = GeoArray(AWAP, :rainfall; date=DateTime(2019, 10, 19))
        @test crs(A) == EPSG(4326)
        st = GeoStack(AWAP; date=DateTime(2019, 10, 19), resize=crop)
        @test crs(st) == EPSG(4326)
        dates = DateTime(2019, 09, 19), DateTime(2019, 11, 19)
        s = GeoSeries(AWAP; date=dates, resize=crop)
        # s = GeoSeries(AWAP; date=dates, resize=resample, crs=EPSG(4326)) TODO: all the same
        # s = GeoSeries(AWAP; date=dates, resize=extend) TODO: this is slow !!!
        @test crs(s[1][:solar]) == EPSG(4326)
        @test A isa GeoArray
        @test st isa GeoStack
        @test s isa GeoSeries
    end
end
