using Rasters, Test, Dates, DimensionalData
using Rasters.LookupArrays, Rasters.Dimensions
using Rasters: FileArray, ASCIIfile

const DD = DimensionalData

ascpath = "britanny.asc"

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
        @test length(val(dims(dims(ascarray), X))) == 401
        @test ndims(ascarray) == 2
        @test dims(ascarray) isa Tuple{<:X,<:Y}
        @test refdims(ascarray) == ()
        @test bounds(ascarray) == ((-4.591421949457, -3.3483663541919997), (47.832347916354, 48.665462015892))
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
        @test ascarray[X(12), Y(56)] == 5340.0
        @test ascarray[35, 42] == 7482.0
        @test ascarray[Y(At(20.0; atol=1e10)), X(At(20; atol=1e10))] == 5500.0
        @test ascarray[Y(Contains(48.5)), X(Contains(-4.2))] == 6624.0
    end

end