using Rasters, Test, GRIBDatasets
using Rasters: FileArray, GRIBfile
using Rasters.LookupArrays, Rasters.Dimensions
using Statistics
using Dates

test_files = [
    "era5-levels-members.grib", # OK
    "fields_with_missing_values.grib", # OK
    "forecast_monthly_ukmo.grib", # Failing because because there are multiple steps
    "lambert_grid.grib", # longitude and latitude are 2D variables. Don't know what to do in this case
    "multi_param_on_multi_dims.grib", # Passing, but we loose the information on number and time. Should be fixed
    "regular_gg_ml.grib", # OK
    "regular_gg_ml_g2.grib", # OK
    "regular_gg_pl.grib", # OK
    "regular_gg_sfc.grib", # OK
    "regular_ll_msl.grib", # OK
    "regular_ll_sfc.grib", # OK
    "regular_ll_wrong_increment.grib", # OK
    "scanning_mode_64.grib", # OK
    "t_analysis_and_fc_0.grib", # OK
]

gribexamples_dir = abspath(joinpath(dirname(pathof(GRIBDatasets)), "..", "test", "sample-data"))

era5 = joinpath(gribexamples_dir, "era5-levels-members.grib")
gribds = GRIBDataset(era5)

@testset "Raster" begin
    @time gribarray = Raster(era5)
    @time lazyarray = Raster(era5; lazy=true);
    @time eagerarray = Raster(era5; lazy=false);
    @time ds = GRIBDataset(era5)

    @testset "lazyness" begin
        @time read(Raster(era5));
        # Eager is the default
        @test parent(gribarray) isa Array
        @test parent(lazyarray) isa FileArray
        @test parent(eagerarray) isa Array
    end

    # @testset "open" begin
    #     @test all(open(A -> A[Y=1], ncarray) .=== ncarray[:, :, :, :, 1])
    # end

    @testset "read" begin
        @time A = read(gribarray);
        @test A isa Raster
        @test parent(A) isa Array
        A2 = zero(A)
        @time read!(gribarray, A2);
        A3 = zero(A)
        @time read!(gribarray, A3)
        @test all(A .=== A2) 
        @test all(A .=== A3)
    end

    @testset "stack, compare to GRIBDataset" begin
        stack = RasterStack(era5)
        ds = GRIBDataset(era5)

        diff = stack[:z][:,:,1,1,1] - ds["z"][:,:,1,1,1]

        @test all(diff .== 0.)
    end

    @testset "array properties" begin
        dsvar = ds["z"]
        @test size(gribarray) == size(dsvar)
        @test gribarray isa Raster
        @test index(gribarray, Ti) == DateTime(2017, 1, 1):Hour(12):DateTime(2017, 1, 2, 12)
        @test index(gribarray, Y) == 90.0:-3.0:-90.0
        @test index(gribarray, X) == 0.0:3.0:357.0
    end

    @testset "dimensions" begin
        @test ndims(gribarray) == 5
        @test length.(dims(gribarray)) == (120, 61, 10, 4, 2)
        @test dims(gribarray) isa Tuple{<:X,<:Y,<:Dim{:number},<:Ti,<:Z}
        @test refdims(gribarray) == ()
        @test typeof(lookup(gribarray)) <: Tuple{<:Mapped,<:Mapped,<:Sampled,<:Sampled,<:Sampled}
        @test bounds(gribarray) == ((0.0, 357.0), (-90.0, 90.0), (0, 9), (DateTime("2017-01-01T00:00:00"), DateTime("2017-01-02T12:00:00")), (500, 850))
    end


    @testset "other fields" begin
        @test ismissing(missingval(gribarray))
        @test metadata(gribarray) isa Metadata{GRIBfile,Dict{String,Any}}
    end

    @testset "indexing" begin
        @test gribarray[Ti(1)] isa Raster{<:Any,4}
        @test gribarray[Y(1), Ti(1)] isa Raster{<:Any,3}
        @test gribarray[X(1), Ti(1)] isa Raster{<:Any,3}
        @test gribarray[X = 2:4, Y = 10:15, Ti = 1, Z = 2, number = 2:4] isa Raster{<:Any,3}
        @test gribarray[X(30), Y(30), Ti(1), Z(1), number = 2] isa Float64
    end

    @testset "methods" begin 
        @testset "mean" begin
            @test all(mean(gribarray; dims=Y) .=== mean(parent(gribarray); dims=2))
        end
        @testset "trim, crop, extend" begin
            a = read(gribarray)
            a[X(1:20)] .= missingval(a)
            trimmed = trim(a)
            @test size(trimmed) == (100, 61, 10, 4, 2)
            cropped = crop(a; to=trimmed)
            @test size(cropped) == (100, 61, 10, 4, 2)
            @test all(collect(cropped .=== trimmed))
            extended = extend(cropped; to=a)
            @test all(collect(extended .=== a))
        end
        @testset "mask and mask!" begin
            msk = read(gribarray)
            msk[X(1:100), Y([1, 5, 45])] .= missingval(msk)
            @test !all(gribarray[X(1:100)] .=== missingval(msk))
            masked = mask(gribarray; with=msk)
            @test all(masked[X(1:100), Y([1, 5, 45])] .=== missingval(msk))
        end
        @testset "slice" begin
            @test_throws ArgumentError Rasters.slice(gribarray, Band)
            ser = Rasters.slice(gribarray, Ti) 
            @test ser isa RasterSeries
            @test size(ser) == (4,)
            @test index(ser, Ti) == DateTime(2017, 1, 1):Hour(12):DateTime(2017, 1, 2, 12)
            @test bounds(ser) == ((DateTime(2017, 1, 1), DateTime(2017, 1, 2, 12)),)
            A = ser[1]
            @test index(A, Y) == 90.0:-3.0:-90.0
            @test index(A, X) == 0.0:3.0:357.0
            @test bounds(A) == ((0.0, 357.0), (-90.0, 90.0), (0, 9), (500, 850))
        end
    end

    # @testset "indexing with reverse lat" begin
    #     if !haskey(ENV, "CI") # CI downloads fail. But run locally
    #         ncrevlat = maybedownload("ftp://ftp.cdc.noaa.gov/Datasets/noaa.ersst.v5/sst.mon.ltm.1981-2010.nc")
    #         ncrevlatarray = Raster(ncrevlat; key=:sst)
    #         @test order(dims(ncrevlatarray, Y)) == ReverseOrdered()
    #         @test ncrevlatarray[Y(At(40)), X(At(100)), Ti(1)] === missing
    #         @test ncrevlatarray[Y(At(-40)), X(At(100)), Ti(1)] === ncrevlatarray[51, 65, 1] == 14.5916605f0
    #         @test val(span(ncrevlatarray, Ti)) == Month(1)
    #         @test val(span(ncrevlatarray, Ti)) isa Month # Not CompoundPeriod
    #     end
    # end

    @testset "selectors" begin
        a = gribarray[X(At(21.0)), Y(Between(50, 52)), Ti(Near(DateTime(2002, 12)))]
        @test bounds(a) == ((51.0, 51.0), (0, 9), (500, 850))
        # x = gribarray[X(Near(150)), Y(Near(30)), Ti(1), number=1, Z(1)]
        # @test x isa Float64
        # lookup(gribarray)
        # dimz = X(Between(-0.0, 360)), Y(Between(-90, 90)), 
        #        Ti(Between(DateTime360Day(2001, 1, 1), DateTime360Day(2003, 01, 02)))
        # @test size(gribarray[dimz...]) == (180, 170, 24)
        # @test index(gribarray[dimz...]) == index(gribarray)
        # nca = gribarray[Y(Between(-80, -25)), X(Between(-0.0, 180.0)), Ti(Contains(DateTime360Day(2002, 02, 20)))]
        # @test size(nca) == (90, 55)
        # @test index(nca, Y) == index(gribarray[1:90, 1:55, 2], Y)
        # @test all(nca .=== gribarray[1:90, 1:55, 14])
    end

    # @testset "conversion to Raster" begin
    #     geoA = ncarray[X(1:50), Y(20:20), Ti(1)]
    #     @test size(geoA) == (50, 1)
    #     @test eltype(geoA) <: Union{Missing,Float32}
    #     @test geoA isa Raster{Union{Missing,Float32},2}
    #     @test dims(geoA) isa Tuple{<:X,<:Y}
    #     @test refdims(geoA) isa Tuple{<:Ti}
    #     @test metadata(geoA) == metadata(ncarray)
    #     @test ismissing(missingval(geoA))
    #     @test name(geoA) == :tos
    # end

    # @testset "show" begin
    #     sh = sprint(show, MIME("text/plain"), ncarray)
    #     # Test but don't lock this down too much
    #     @test occursin("Raster", sh)
    #     @test occursin("Y", sh)
    #     @test occursin("X", sh)
    #     @test occursin("Time", sh)
    # end

    # @testset "plot" begin
    #     ncarray[Ti(1:3:12)] |> plot
    #     ncarray[Ti(1)] |> plot
    #     ncarray[Y(100), Ti(1)] |> plot
    # end

end