using Rasters, Test, GRIBDatasets
using Rasters: FileArray
using Rasters.LookupArrays, Rasters.Dimensions
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

gribfiles = joinpath.(gribexamples_dir, test_files)
era5 = joinpath(gribexamples_dir, "era5-levels-members.grib")

# for test_file in gribfiles
#     print("Testing: ", test_file)
#     stack = try
#         RasterStack(test_file)
#         println(lpad("OK", 5))
#     catch
#         println(lpad("NOK", 5))
#     end
# end

@testset "Raster" begin
    @time ncarray = Raster(era5)
    @time lazyarray = Raster(era5; lazy=true);
    @time eagerarray = Raster(era5; lazy=false);
    @time ds = GRIBDataset(era5)

    @testset "lazyness" begin
        @time read(Raster(era5));
        # Eager is the default
        @test parent(ncarray) isa Array
        @test parent(lazyarray) isa FileArray
        @test parent(eagerarray) isa Array
    end

    # @testset "open" begin
    #     @test all(open(A -> A[Y=1], ncarray) .=== ncarray[:, :, :, :, 1])
    # end

    @testset "read" begin
        @time A = read(ncarray);
        @test A isa Raster
        @test parent(A) isa Array
        A2 = zero(A)
        @time read!(ncarray, A2);
        A3 = zero(A)
        @time read!(ncarray, A3)
        @test all(A .=== A2) 
        @test all(A .=== A3)
    end

    @testset "stack" begin
        stack = RasterStack(era5)
        ds = GRIBDataset(era5)

        diff = stack[:z][:,:,1,1,1] - ds["z"][:,:,1,1,1]

        @test all(diff .== 0.)
    end

    # @testset "array properties" begin
    #     dsvar = ds["z"]
    #     @test size(ncarray) == size(dsvar)
    #     @test ncarray isa Raster
    #     @test index(ncarray, Ti) == DateTime(2017, 1, 1):Hour(12):DateTime(2017, 1, 2, 12)
    #     @test index(ncarray, Y) == -79.5:89.5
    #     @test index(ncarray, X) == 1.0:2:359.0
    #     @test bounds(ncarray) == (
    #         (0.0, 360.0), 
    #         (-80.0, 90.0), 
    #         (DateTime360Day(2001, 1, 1), DateTime360Day(2003, 1, 1)),
    #     )
    # end

    @testset "dimensions" begin
        @test ndims(ncarray) == 5
        @test length.(dims(ncarray)) == (120, 61, 10, 4, 2)
        @test dims(ncarray) isa Tuple{<:X,<:Y,<:Dim{:number},<:Ti,<:Z}
        @test refdims(ncarray) == ()
        # @test val.(span(ncarray)) == 
        #     (vcat((0.0:2.0:358.0)', (2.0:2.0:360.0)'),
        #      vcat((-80.0:89.0)', (-79.0:90.0)'),
        #      vcat(permutedims(DateTime360Day(2001, 1, 1):Month(1):DateTime360Day(2002, 12, 1)), 
        #           permutedims(DateTime360Day(2001, 2, 1):Month(1):DateTime360Day(2003, 1, 1)))
        #     )
        # this test not passing 
        @test typeof(lookup(ncarray)) <: Tuple{<:Mapped,<:Mapped,<:Sampled,<:Sampled,<:Sampled}
        @test bounds(ncarray) == ((0.0, 357.0), (-90.0, 90.0), (0, 9), (DateTime("2017-01-01T00:00:00"), DateTime("2017-01-02T12:00:00")), (500, 850))
    end
    # tempfile = tempname() * ".nc"


    # @testset "other fields" begin
    #     @test ismissing(missingval(ncarray))
    #     @test metadata(ncarray) isa Metadata{NCDfile}
    #     @test name(ncarray) == :tos
    # end

    @testset "indexing" begin
        @test ncarray[Ti(1)] isa Raster{<:Any,4}
        @test ncarray[Y(1), Ti(1)] isa Raster{<:Any,3}
        @test ncarray[X(1), Ti(1)] isa Raster{<:Any,3}
        # @test ncarray[X(1), Y(1), Ti(1)] isa Missing
        @test ncarray[X(30), Y(30), Ti(1), Z(1), number = 2] isa Float64
        # Russia
        # @test ncarray[X(50), Y(100), Ti(1)] isa Missing
        # Alaska
        # @test ncarray[Y(Near(64.2008)), X(Near(149.4937)), Ti(1)] isa Missing
        # @test ncarray[Ti(2), X(At(59.0)), Y(At(-50.5))] == ncarray[30, 30, 2] === 278.47168f0
    end

    # @testset "methods" begin 
    #     @testset "mean" begin
    #         @test all(mean(ncarray; dims=Y) .=== mean(parent(ncarray); dims=2))
    #     end
    #     @testset "trim, crop, extend" begin
    #         a = read(ncarray)
    #         a[X(1:20)] .= missingval(a)
    #         trimmed = trim(a)
    #         @test size(trimmed) == (160, 169, 24)
    #         cropped = crop(a; to=trimmed)
    #         @test size(cropped) == (160, 169, 24)
    #         @test all(collect(cropped .=== trimmed))
    #         extended = extend(cropped; to=a)
    #         @test all(collect(extended .=== a))
    #     end
    #     @testset "mask and mask!" begin
    #         msk = read(ncarray)
    #         msk[X(1:100), Y([1, 5, 95])] .= missingval(msk)
    #         @test !all(ncarray[X(1:100)] .=== missingval(msk))
    #         masked = mask(ncarray; with=msk)
    #         @test all(masked[X(1:100), Y([1, 5, 95])] .=== missingval(msk))
    #         tempfile = tempname() * ".nc"
    #         cp(gribfile, tempfile)
    #         @test !all(Raster(tempfile)[X(1:100), Y([1, 5, 95])] .=== missing)
    #         open(Raster(tempfile; lazy=true); write=true) do A
    #             mask!(A; with=msk, missingval=missing)
    #             # TODO: replace the CFVariable with a FileArray{NCDfile} so this is not required
    #             nothing
    #         end
    #         @test all(Raster(tempfile)[X(1:100), Y([1, 5, 95])] .=== missing)
    #         rm(tempfile)
    #     end
    #     @testset "mosaic" begin
    #         @time ncarray = Raster(gribfile)
    #         A1 = ncarray[X(1:80), Y(1:100)]
    #         A2 = ncarray[X(50:150), Y(90:150)]
    #         tempfile = tempname() * ".nc"
    #         Afile = mosaic(first, read(A1), read(A2); missingval=missing, atol=1e-7, filename=tempfile)
    #         Amem = mosaic(first, A1, A2; missingval=missing, atol=1e-7)
    #         Atest = ncarray[X(1:150), Y(1:150)]
    #         Atest[X(1:49), Y(101:150)] .= missing
    #         Atest[X(81:150), Y(1:89)] .= missing
    #         @test all(Atest .=== Afile .=== Amem)
    #     end
    #     @testset "slice" begin
    #         @test_throws ArgumentError Rasters.slice(ncarray, Z)
    #         ser = Rasters.slice(ncarray, Ti) 
    #         @test ser isa RasterSeries
    #         @test size(ser) == (24,)
    #         @test index(ser, Ti) == DateTime360Day(2001, 1, 16):Month(1):DateTime360Day(2002, 12, 16)
    #         @test bounds(ser) == ((DateTime360Day(2001, 1, 1), DateTime360Day(2003, 1, 1)),)
    #         A = ser[1]
    #         @test index(A, Y) == -79.5:89.5
    #         @test index(A, X) == 1.0:2:359.0
    #         @test bounds(A) == ((0.0, 360.0), (-80.0, 90.0))
    #     end
    # end

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

    # @testset "selectors" begin
    #     a = ncarray[X(At(21.0)), Y(Between(50, 52)), Ti(Near(DateTime360Day(2002, 12)))]
    #     @test bounds(a) == ((50.0, 52.0),)
    #     x = ncarray[X(Near(150)), Y(Near(30)), Ti(1)]
    #     size(ncarray)
    #     @test x isa Float32
    #     lookup(ncarray)
    #     dimz = X(Between(-0.0, 360)), Y(Between(-90, 90)), 
    #            Ti(Between(DateTime360Day(2001, 1, 1), DateTime360Day(2003, 01, 02)))
    #     @test size(ncarray[dimz...]) == (180, 170, 24)
    #     @test index(ncarray[dimz...]) == index(ncarray)
    #     nca = ncarray[Y(Between(-80, -25)), X(Between(-0.0, 180.0)), Ti(Contains(DateTime360Day(2002, 02, 20)))]
    #     @test size(nca) == (90, 55)
    #     @test index(nca, Y) == index(ncarray[1:90, 1:55, 2], Y)
    #     @test all(nca .=== ncarray[1:90, 1:55, 14])
    # end

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