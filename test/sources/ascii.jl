using Rasters, Test, Dates, DimensionalData, RasterDataSources, Statistics, Plots
using Rasters.LookupArrays, Rasters.Dimensions
using Rasters: FileArray, ASCIIfile

const DD = DimensionalData

ascpath = getraster(MOD13Q1, :NDVI; RasterDataSources.crozon...)

@testset "ASCII array" begin
    @time ascarray = Raster(ascpath)
    @time lazyarray = Raster(ascpath; lazy=true)
    @time eagerarray = Raster(ascpath; lazy=false)
    @test_throws ArgumentError Raster("notafile.asc")

    @testset "lazyness" begin
        @test parent(ascarray) isa Array
        @test parent(lazyarray) isa FileArray
        @test parent(eagerarray) isa Array
    end

    @testset "open" begin
        @test all(open(A -> A[Y=1], ascarray) .=== ascarray[:, 1, :])
        tempfile = tempname()
        cp(ascpath, tempfile * ".asc")
        ascwritearray = Raster(tempfile * ".asc"; lazy=true)
        open(ascwritearray; write=true) do A
            A .*= 2
        end
        @test Raster(tempfile * ".asc") == ascarray .* 2
    end

    @testset "read" begin
        A = read(ascarray)
        @test A isa Raster
        @test parent(A) isa Array
        A2 = zero(A)
        @time read!(ascarray, A2);
        A3 = zero(A)
        @time read!(ascpath, A3);
        @test A == A2 == A3
    end
        
    @testset "array properties" begin
        @test ascarray isa Raster{Float64,2}
    end

    @testset "dimensions" begin
        @test length(val(dims(dims(ascarray), X))) == 9
        @test ndims(ascarray) == 2
        @test dims(ascarray) isa Tuple{<:X,<:Y}
        @test refdims(ascarray) == ()
        @test bounds(ascarray) == ((-4.513046935220226, -4.484930533174481), (48.231304357467636, 48.24592039430238))
    end

    @testset "fields" begin
        @test missingval(ascarray) == -9999.0
        @test metadata(ascarray) isa Metadata{ASCIIfile}
        @test name(ascarray) == Symbol("")
        @test label(ascarray) == ""
        @test units(ascarray) === nothing
        custom = Raster(ascpath; name = :ascii, mappedcrs = EPSG(4326))
        @test name(custom) == :ascii
        @test label(custom) == "ascii"
        @test mappedcrs(dims(custom, Y)) == EPSG(4326)
        @test mappedcrs(dims(custom, X)) == EPSG(4326)
        @test mappedcrs(custom) == EPSG(4326)
        @test crs(dims(custom, Y)) == EPSG(4326)
        @test crs(dims(custom, X)) == EPSG(4326)
        @test crs(custom) == EPSG(4326)
    end

    @testset "getindex" begin
        @test ascarray[X(1)] isa Raster{Float64, 1}
        @test ascarray[Y(1)] isa Raster{Float64, 1}
        @test ascarray[X(8), Y(7)] == 2425.0
        @test ascarray[6, 3] == 7432.0
        @test ascarray[Y(At(1; atol=1e10)), X(At(2; atol=1e10))] == -3000.0
        @test ascarray[Y(Contains(48.24)), X(Contains(-4.49))] == 5592.0
    end

    @testset "methods" begin
        @test mean(ascarray) == 5736.0
        @test mean(ascarray; dims=Y) == mean(parent(ascarray); dims=2)

        @testset "trim, crop, extend" begin
            a = read(ascarray)
            a[X(1:2)] .= missingval(a)
            trimmed = trim(a)
            @test size(trimmed) == (7, 9)
            cropped = crop(a; to = trimmed)
            @test size(cropped) == (7, 9)
            # @test all(collect(cropped .== trimmed))
            # extended = extend(cropped; to = b)
            # @test all(collect(extended .== b))
        end
    end

end