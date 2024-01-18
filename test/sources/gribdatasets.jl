using Rasters, Test, GRIBDatasets
using Rasters: FileArray, CDMsource
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

@testset "Raster" begin
    @time gribarray = Raster(era5)
    @time lazyarray = Raster(era5; lazy=true);
    @time lazystack= RasterStack(era5; lazy=true);
    @time eagerstack = RasterStack(era5; lazy=false);
    @time ds = GRIBDataset(era5)

    @testset "lazyness" begin
        @time read(Raster(era5));
        @test parent(gribarray) isa Array
        @test parent(lazyarray) isa FileArray
    end

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
        stack = RasterStack(era5; lazy = true)
        ds = GRIBDataset(era5)

        diff = stack[:z][:,:,1,1,1] - ds["z"][:,:,1,1,1]

        @test all(diff .== 0.)
    end

    @testset "eager stack" begin
        t = eagerstack[:t]
        @test t[:,:,2,3,1] isa AbstractMatrix
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
        @test length.(dims(gribarray)) == (120, 61, 2, 10, 4)
        @test dims(gribarray) isa Tuple{<:X,<:Y,<:Z,<:Dim{:number},<:Ti}
        @test refdims(gribarray) == ()
        @test typeof(lookup(gribarray)) <: Tuple{<:Mapped,<:Mapped,<:Sampled,<:Sampled,<:Sampled}
        @test Rasters.bounds(gribarray) == ((0.0, 357.0), (-90.0, 90.0), (500, 850), (0, 9), (DateTime("2017-01-01T00:00:00"), DateTime("2017-01-02T12:00:00")))
    end

    @testset "cf attributes" begin
        z = lazystack[:z]
        @test metadata(z)["standard_name"] == "geopotential"

        @test metadata(lazystack)["Conventions"] == "CF-1.7"
        x = dims(lazystack, :X)
        @test metadata(x)["standard_name"] == "longitude"
    end

    @testset "other fields" begin
        @test ismissing(missingval(gribarray))
        @test metadata(gribarray) isa Metadata{CDMsource,Dict{String,Any}}
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
            @test size(trimmed) == (100, 61, 2, 10, 4)
            cropped = crop(a; to=trimmed)
            @test size(cropped) == (100, 61, 2, 10, 4)
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
            @test_throws DimensionMismatch Rasters.slice(gribarray, Band)
            ser = Rasters.slice(gribarray, Ti) 
            @test ser isa RasterSeries
            @test size(ser) == (4,)
            @test index(ser, Ti) == DateTime(2017, 1, 1):Hour(12):DateTime(2017, 1, 2, 12)
            @test Rasters.bounds(ser) == ((DateTime(2017, 1, 1), DateTime(2017, 1, 2, 12)),)
            A = ser[1]
            @test index(A, Y) == 90.0:-3.0:-90.0
            @test index(A, X) == 0.0:3.0:357.0
        end
    end

    @testset "selectors" begin
        a = gribarray[X(At(21.0)), Y(Between(50, 52)), Ti(Near(DateTime(2002, 12)))]
        @test Rasters.bounds(a) == ((51.0, 51.0), (500, 850), (0, 9))
    end
end