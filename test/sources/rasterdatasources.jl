using Rasters, RasterDataSources, Test, Dates, ArchGDAL, NCDatasets, Proj

# Too big to test on CI
# if !haskey(ENV, "CI")
    # @testset "load WorldClim Weather" begin
    #     # Weather time-series
    #     ser = RasterSeries(WorldClim{Weather}, (:prec,); 
    #         date=(Date(2001), Date(2002)), missingval=NaN32
    #     )
    #     @test all(ser[At(Date(2001, 1))].prec .===
    #         Raster(WorldClim{Weather}, :prec; date=DateTime(2001), missingval=NaN32))
    # end
# end

@testset "load WorldClim Climate" begin
    # Weather time-series
    ser = RasterSeries(WorldClim{Climate}, :prec; res="10m", month=Jan:March, mappedcrs=EPSG(4326), raw=true)
    # Select Australia, using regular lat/lon selectors
    A = ser[month=Jan]
    @test A isa Raster
    A[Y(Between(-10, -45)), X(Between(110, 160))]
    st = RasterStack(WorldClim{Climate}, (:prec, :tmax); month=1, raw=true)
    @test st[:prec] == A
    @test missingval(st) == (prec=-32768, tmax=-3.4f38)
    @test st isa RasterStack{(:prec,:tmax),@NamedTuple{prec::Int16,tmax::Float32},2}
end

@testset "load WorldClim BioClim" begin
    A = Raster(WorldClim{BioClim}, :Bio_1; mappedcrs=EPSG(4326))
    A[Y(Between(-10, -45)), X(Between(110, 160))]
    @test A isa Raster
    @test missingval(A) === missing
    st = RasterStack(WorldClim{BioClim}, (1, 2))
    @test all(st.bio1 .=== A)
    @test st isa RasterStack
    @test A isa Raster
    # Future Bioclim works
    st = RasterStack(WorldClim{Future{BioClim, CMIP6, GFDL_ESM4, SSP370}}, (1, 2); 
        date = Date(2050), res = "10m",
        lazy=true, 
        missingval=Inf32, 
        crs=nothing, 
        mappedcrs=EPSG(4326),
    )
    @test missingval(st) === Inf32
    @test missingval(st.bio1) === Inf32
    ra = Raster(WorldClim{Future{BioClim, CMIP6, GFDL_ESM4, SSP370}}, 2; 
        date = Date(2050), res = "10m"
    )
    @test Rasters.name(ra) == :bio2    
end

@testset "load CHELSA BioClim" begin
    A = Raster(CHELSA{BioClim}, 1; lazy=true)
    @test Rasters.name(A) == :bio1
    st = RasterStack(CHELSA{BioClim}, (:bio1, :BIO2); lazy=true)
    @test keys(st) == (:bio1, :bio2)
    @test A isa Raster{Float64,2}
    @test st isa RasterStack
    @test st.bio2 isa Raster{Float64,2}

    A = Raster(CHELSA{BioClim}, 1; lazy=true, raw=true)
    st = RasterStack(CHELSA{BioClim}, (:bio1, :BIO2); lazy=true, raw=true)
    @test A isa Raster{UInt16,2}
    @test st isa RasterStack{(:bio1, :bio2),@NamedTuple{bio1::UInt16, bio2::UInt16}}
    @test st.bio2 isa Raster{UInt16,2}

    st = RasterStack(CHELSA{BioClim}, (1, 2); lazy=true)

    @test missingval(st) === missingval(st.bio1) === nothing
    @test metadata(st) == Rasters.NoMetadata()
end

@testset "load EarthEnv HabitatHeterogeneity" begin
    A = Raster(EarthEnv{HabitatHeterogeneity}, :cv; mappedcrs=EPSG(4326))
    A[Y(Between(-10, -45)), X(Between(110, 160))] 
    st = RasterStack(EarthEnv{HabitatHeterogeneity}, (:cv, :evenness))
    @test A isa Raster
    @test st isa RasterStack
    @test st[:evenness] isa Raster
end

@testset "load EarthEnv LandCover" begin
    A = Raster(EarthEnv{LandCover}, 2; mappedcrs=EPSG(4326))
    @test Rasters.name(A) == :evergreen_broadleaf_trees
    A[Y(Between(-10, -45)), X(Between(110, 160))]
    @test A isa Raster
    st = RasterStack(EarthEnv{LandCover}, (:evergreen_broadleaf_trees, :deciduous_broadleaf_trees); mappedcrs=EPSG(4326))
    @test keys(st) == (:evergreen_broadleaf_trees, :deciduous_broadleaf_trees)
    @test st isa RasterStack
end

@testset "load ALWB" begin
    A = Raster(ALWB{Deciles,Day}, :rain_day; date=DateTime(2019, 10, 19), lazy=true)
    @test crs(A) == EPSG(4326)
    A = Raster(ALWB{Values,Day}, :ss_pct; date=DateTime(2019, 10, 19), lazy=true)
    @test crs(A) == EPSG(4326)
    st = RasterStack(ALWB{Values,Day}, (:s0_pct, :ss_pct); date=DateTime(2019, 10, 19), lazy=true)
    @test crs(st) == EPSG(4326)
    @test crs(st[:s0_pct]) == EPSG(4326)
    dates = DateTime(2019, 10, 19), DateTime(2021, 11, 20)
    s = RasterSeries(ALWB{Values,Day}, (:s0_pct, :ss_pct); date=dates, lazy=true)
    @test A isa Raster
    @test st isa RasterStack
    @test s isa RasterSeries
end

# Obscure .Z format may not work on windows
if Sys.islinux()
    @testset "load AWAP" begin
        A = Raster(AWAP, :rainfall; date=DateTime(2019, 10, 19), lazy=true)
        @test crs(A) == EPSG(4326)
        # ALWB :solar has a broken index - the size is different and the
        # points don't exactly match the other layers.
        # Need to work out how to best resolve this kind of problem so that we can
        # still use the layers in stacks.
        layers = (:rainfall, :vprpress09, :vprpress15, :tmin, :tmax)
        st = RasterStack(AWAP, layers; date=DateTime(2019, 10, 19), resize=crop, lazy=true)
        @test crs(st) == EPSG(4326)
        dates = DateTime(2019, 09, 19), DateTime(2019, 11, 19)
        s = RasterSeries(AWAP, layers; date=dates, resize=crop, lazy=true)
        # test date as an Array
        s2 = RasterSeries(AWAP, layers; date=[dates...], resize=crop)
        # s = RasterSeries(AWAP; date=dates, resize=resample, crs=EPSG(4326)) TODO: all the same
        # s = RasterSeries(AWAP; date=dates, resize=extend) TODO: this is slow !!!
        @test crs(s[1][:rainfall]) == EPSG(4326)
        @test A isa Raster
        @test st isa RasterStack
        @test s isa RasterSeries
        @test s2 isa RasterSeries
        @test length(s2) == 2 # date is an array: don't take intermediate steps
    end
end
