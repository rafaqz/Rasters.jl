using GeoData, Test, Dates, Statistics
using GeoData: Start, Center, End,
      formatdims, dims, aggregate, upsample, downsample


@testset "upsample" begin
    @test upsample(1, 2) == 1
    @test upsample(2, 2) == 3
    @test upsample(3, 2) == 5
    @test upsample(1, 3) == 1
    @test upsample(2, 3) == 4
    @test upsample(3, 3) == 7
end

@testset "downsample" begin
    @test downsample(1, 2) == 1
    @test downsample(3, 2) == 2
    @test downsample(5, 2) == 3
    @test downsample(1, 3) == 1
    @test downsample(4, 3) == 2
    @test downsample(7, 3) == 3
end

data1 = [ 1  2  3  4  5  6 -1
          7  8  9 10 11 12 -1
         13 14 15 16 17 18 -1]
data2 = 2 * data1
data3 = 3 * data1
data4 = 4 * data1
dimz = Lon([30, 40, 50]), Lat(LinRange(-10, 24, 7))
array1 = GeoArray(data1, dimz)
array2 = GeoArray(data2, dimz)
array1a = GeoArray(data3, dimz)
array2a = GeoArray(data4, dimz)
stack1 = GeoStack(array1, array2; keys=(:array1, :array2))
stack2 = GeoStack(array1a, array2a; keys=(:array1, :array2))
dates = [DateTime(2017), DateTime(2018)]
series = GeoSeries([stack1, stack2], (Ti(dates),));

@testset "Aggregate at a locus" begin
    @testset "single scale single locus" begin
        scale = 3
        @test aggregate(array1, Start(), scale) == [1 4]
        @test aggregate(array1, Center(), scale) == [8 11]
        @test aggregate(array1, End(), scale) == [15 18]
        @test aggregate(stack1, Start(), scale)[:array2] == [2 8]
        @test aggregate(stack1, Center(), scale)[:array2] == [16 22]
        @test map(x -> aggregate(x, Start(), scale), series)[2][:array2] == [4 16]
        @test typeof(map(x -> aggregate(x, Start(), scale), series)) <: GeoSeries
        A = aggregate(array1, Start(), scale)
        @test length.(dims(A)) == size(A)
    end
    @testset "mixed scales" begin
        scale = (3, 2)
        @test aggregate(array1, Start(), scale) == [1 3 5]
        scale = (Lat(2), Lon(3))
        @test aggregate(array1, Center(), scale) == [8 10 12]
        @test aggregate(array1, End(), scale) == [14 16 18]
        A = aggregate(array1, Start(), scale)
        @test length.(dims(A)) == size(A)
    end
    @testset "mixed locus" begin
        scale = 3
        @test aggregate(array1, (End(), Start()), 3) == [13 16]
        @test aggregate(array1, (End(), Start()), (3, 2)) == [13 15 17]
        A = aggregate(array1, (End(), Start()), scale)
        @test length.(dims(A)) == size(A)
    end
end

@testset "Aggregate with a function" begin
    @test aggregate(array1, sum, 3) == [72 99]
    @test aggregate(array1, median, 3) == [8 11]
    @test aggregate(array1, sum, (3, 2)) == [45 57 69]
    A = aggregate(array1, sum, (3, 2))
    @test length.(dims(A)) == size(A)
end
