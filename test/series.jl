using Rasters, Test, Dates, DimensionalData
using Rasters.Lookups, Rasters.Dimensions
using Rasters: metadata
include(joinpath(dirname(pathof(Rasters)), "../test/test_utils.jl"))

# RasterSeries from Raster/RasterStack components

data1 = [1 2 3 4
         5 6 7 8]
data2 = 2 * data1
data3 = 3 * data1
data4 = 4 * data1
dimz = X([30, 40]), Y(-10.0:10.0:20.0)
r1 = Raster(data1, dimz)
r2 = Raster(data2, dimz)
r1a = Raster(data3, dimz)
r2a = Raster(data4, dimz)
stack1 = RasterStack(r1, r2; name=(:r1, :r2))
stack2 = RasterStack(r1a, r2a; name=(:r1, :r2))
dates = [DateTime(2017), DateTime(2018)]
ser = RasterSeries([stack1, stack2], Ti(dates))
@test issorted(dates)

@testset "getindex returns the currect types" begin
    @test ser[Ti(1)] isa RasterStack{(:r1, :r2),@NamedTuple{r1::Int, r2::Int},2,<:NamedTuple}
    @test ser[Ti(1)][:r2] isa Raster{Int,2}
    @test ser[Ti(1)][:r2][1, 1] isa Int
end

@testset "map" begin
    @test map(x -> x[1], ser) == Raster([(r1=1, r2=2), (r1=3, r2=4)], dims(ser))
    @test map(x -> x, ser) == ser;
end

@testset "properties" begin
    @test refdims(ser) === ()
    # Should these be real fields? what is the use-case?  
    @test metadata(ser) === NoMetadata()
    @test name(ser) === DimensionalData.NoName()
    @test label(ser) === ""
end

@testset "getindex returns the currect results" begin
    @test ser[Ti(Near(DateTime(2017)))][:r1][X(1), Y(3)] === 3
    @test ser[Ti(At(DateTime(2018)))][:r2][X(2), Y(4)] === 32
    @test ser[Ti(1)][:r1][X(1), Y(2)] == 2
    @test ser[Ti(1)][:r2][X(2), Y(3:4)] == [14, 16] 
end

@testset "getindex is type stable all the way down" begin
    # @inferred ser[Ti(At(DateTime(2017)))][:r1, X(1), Y(2)]
    @inferred ser[Ti(1)][:r1][X(1), Y(2)]
    # @inferred ser[Ti(1)][:r1, X(1), Y(2:4)]
    # @inferred ser[Ti(1)][:r1][X(1), Y(2:4)]
    # @inferred ser[1][:r1, X(1:2), Y(:)]
    # @inferred ser[1][:r1][X(1:2), Y(:)]
end

@testset "setindex!" begin
    ser[1] = ser[2]
    @test typeof(ser[1]) == typeof(ser[2]) == eltype(parent(ser))
    typeof(parent(ser)[1]) == eltype(ser)
    ser[Ti(1)] = ser[Ti(2)]
    @test ser[Ti(1)] == ser[Ti(2)]
end

@testset "rebuild" begin
    @test rebuild(ser, parent(ser)) === ser
    @test rebuild(ser; dims=dims(ser)) === ser
    @test rebuild(ser; dims=(X(),)) !== ser
    @test rebuild(ser; name=nothing) === ser
    @test rebuild(ser; metadata=nothing) === ser
end

@testset "slice, combine" begin
    r1 = Raster(ones(4, 5, 10), (X(), Y(), Ti(10:10:100))) .* reshape(1:10, (1, 1, 10))
    r2 = r1 .* 2
    ser = slice(r1, Ti)
    @test size(ser) == (10,)
    combined = Rasters.combine(ser, Ti())
    @test combined == r1
    @test dims(combined) == dims(r1)
    ser = slice(r1, (X, Ti))
    @test size(ser) == (4, 10)
    ser = slice(r1, (X, Y, Ti))
    combined2 = Rasters.combine(ser, (X, Y, Ti))
    slice(r1, Ti)
    @test combined == r1 == permutedims(combined2, (X, Y, Ti))
    @test dims(combined) == dims(r1) == dims(permutedims(combined2, (X, Y, Ti)))
    stack = RasterStack((r1=r1, r2=r2))
    ser = slice(stack, Ti)
    @test size(ser) == (10,)
    combined = Rasters.combine(ser, Ti)
    dims(first(combined))
    ser = slice(stack, (Y, Ti))
    @test size(ser) == (5, 10,)
    combined = Rasters.combine(ser, (Y, Ti))
end

@testset "show" begin
    # 2d
    ser2 = slice(r1, (X, Y))
    sh = sprint(show, MIME("text/plain"), ser2)
    # Test but don't lock this down too much
    @test occursin("RasterSeries", sh)
    @test occursin("Raster", sh)
    @test occursin("X", sh)
    @test occursin("Y", sh)
    # 1d
    ser1 = slice(r1, X)
    sh = sprint(show, MIME("text/plain"), ser1)
    # Test but don't lock this down too much
    @test occursin("RasterSeries", sh)
    @test occursin("Raster", sh)
    @test occursin("X", sh)
end

@testset "duplicate_first & lazy" begin
    temporary_random_rasters(3, (10,10,3), UInt8) do filenames
        times = Ti(DateTime.(sort(rand(UInt16, length(filenames)))))
        series = RasterSeries(filenames, times; duplicate_first=true, lazy=true)
        @test all(Rasters.filename.(series) .== filenames)
        first_dims = dims(first(series))
        @test all(dims(r) == first_dims for r in series)
        @test Rasters.isdisk(series)
        @test !Rasters.isdisk(read(series))
        @test Rasters.isdisk(Rasters.combine(series))
        end
end
