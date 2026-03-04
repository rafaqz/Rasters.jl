using Rasters, RasterDataSources, Test, Dates, ArchGDAL, NCDatasets, Proj, Missings
using Rasters.DiskArrays

# Skip slow/large downloads on CI unless explicitly enabled
const RUN_LARGE_TESTS = get(ENV, "RASTERS_RUN_LARGE_TESTS", "false") == "true"

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

# ALWB data source URLs have changed and are no longer available at old URLs
# Skipping tests until RasterDataSources is updated with new URLs
# @testset "load ALWB" begin
#     A = Raster(ALWB{Deciles,Day}, :rain_day; date=DateTime(2019, 10, 19), lazy=true)
#     @test crs(A) == EPSG(4326)
#     A = Raster(ALWB{Values,Day}, :ss_pct; date=DateTime(2019, 10, 19), lazy=true)
#     @test crs(A) == EPSG(4326)
#     st = RasterStack(ALWB{Values,Day}, (:s0_pct, :ss_pct); date=DateTime(2019, 10, 19), lazy=true)
#     @test crs(st) == EPSG(4326)
#     @test crs(st[:s0_pct]) == EPSG(4326)
#     dates = DateTime(2019, 10, 19), DateTime(2021, 11, 20)
#     s = RasterSeries(ALWB{Values,Day}, (:s0_pct, :ss_pct); date=dates, lazy=true)
#     @test A isa Raster
#     @test st isa RasterStack
#     @test s isa RasterSeries
# end

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

@testset "load SRTM" begin
    # Single tile - small area in Switzerland
    A = Raster(SRTM; bounds=(7.0, 8.0, 46.0, 47.0), lazy=true)
    @test A isa Raster
    @test crs(A) == EPSG(4326)
    @test ndims(A) == 2
    @test hasdim(A, X)
    @test hasdim(A, Y)
    # Test that data is accessible (SRTM uses Int16 with missing values)
    @test Missings.nonmissingtype(eltype(A)) <: Integer
    # Test lazy vs eager loading
    A_eager = Raster(SRTM; bounds=(7.0, 8.0, 46.0, 47.0), lazy=false)
    @test !(parent(A_eager) isa DiskArrays.AbstractDiskArray)

    # Multi-tile test (2x2 tiles) - larger area
    if RUN_LARGE_TESTS
        A_multi = Raster(SRTM; bounds=(5.0, 15.0, 45.0, 50.0), lazy=true)
        @test A_multi isa Raster
        @test crs(A_multi) == EPSG(4326)
        # Should be larger than single tile
        @test size(A_multi, X) > size(A, X)
        @test size(A_multi, Y) > size(A, Y)
    end
end

@testset "load TerraClimate" begin
    # TerraClimate returns NetCDF files with monthly time dimension
    A = Raster(TerraClimate{Historical}, :tmax; date=Date(2020), lazy=true)
    @test A isa Raster
    @test ndims(A) == 3  # X, Y, Ti (12 months)
    @test hasdim(A, Ti)
    @test size(A, Ti) == 12

    st = RasterStack(TerraClimate{Historical}, (:tmax, :tmin); date=Date(2020), lazy=true)
    @test st isa RasterStack
    @test keys(st) == (:tmax, :tmin)
end

@testset "load NCEP" begin
    # Test daily pressure data
    A = Raster(NCEP{DailyPressure}, :air; date=Date(2020), lazy=true)
    @test A isa Raster
    @test hasdim(A, Ti)  # Daily data has time dimension

    # Test surface data
    A_surf = Raster(NCEP{SurfaceGauss}, :tmax; date=Date(2020), lazy=true)
    @test A_surf isa Raster

    st = RasterStack(NCEP{DailyPressure}, (:air, :hgt); date=Date(2020), lazy=true)
    @test st isa RasterStack
    @test keys(st) == (:air, :hgt)
end

@testset "load CHELSA Future" begin
    # CHELSA Future CMIP6
    A = Raster(CHELSA{Future{BioClim, CMIP6, GFDL_ESM4, SSP370}}, 1; date=Date(2050), lazy=true)
    @test A isa Raster
    @test Rasters.name(A) == :bio1

    st = RasterStack(CHELSA{Future{BioClim, CMIP6, GFDL_ESM4, SSP370}}, (1, 2); date=Date(2050), lazy=true)
    @test st isa RasterStack
    @test keys(st) == (:bio1, :bio2)

    # CHELSA Future CMIP5 (different date range: 2041-2060 or 2061-2080)
    if RUN_LARGE_TESTS
        A_cmip5 = Raster(CHELSA{Future{BioClim, CMIP5, ACCESS1_0, RCP45}}, 1; date=Date(2050), lazy=true)
        @test A_cmip5 isa Raster
    end
end

# MODIS tests require network access and specific coordinates
# These are slower due to API calls
if RUN_LARGE_TESTS
    @testset "load MODIS" begin
        # Test single date, single layer
        A = Raster(MOD13Q1, :NDVI;
            lat=48.0, lon=-4.0, km_ab=10, km_lr=10,
            date=Date(2020, 6, 1)
        )
        @test A isa Raster
        @test Rasters.name(A) == :NDVI

        # Test stack with multiple layers
        st = RasterStack(MOD13Q1, (:NDVI, :EVI);
            lat=48.0, lon=-4.0, km_ab=10, km_lr=10,
            date=Date(2020, 6, 1)
        )
        @test st isa RasterStack
        @test keys(st) == (:NDVI, :EVI)

        # Test date range returns series
        ser = Raster(MOD13Q1, :NDVI;
            lat=48.0, lon=-4.0, km_ab=10, km_lr=10,
            date=(Date(2020, 1, 1), Date(2020, 3, 1))
        )
        @test ser isa RasterSeries
    end
end
