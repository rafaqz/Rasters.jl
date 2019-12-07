using GeoData, Test, Dates, Statistics
using GeoData: Time, formatdims, dims, aggregate, Start, Center, End

data1 = [ 1  2  3  4  5  6
          7  8  9 10 11 12    
         13 14 15 16 17 18]
data2 = 2 * data1
data3 = 3 * data1
data4 = 4 * data1
dimz = Lon([30, 40, 50]), Lat((-10, 22))
array1 = GeoArray(data1, dimz)
array2 = GeoArray(data2, dimz)
array1a = GeoArray(data3, dimz)
array2a = GeoArray(data4, dimz)
stack1 = GeoStack(array1, array2; keys=(:array1, :array2))
stack2 = GeoStack(array1a, array2a; keys=(:array1, :array2))
dates =[DateTime(2017), DateTime(2018)]
series = GeoSeries([stack1, stack2], (Time(dates),));

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
    end
    @testset "mixed scales" begin
        scale = (3, 2)
        @test aggregate(array1, Start(), scale) == [1 3 5]
        scale = (Lat(2), Lon(3))
        @test aggregate(array1, Center(), scale) == [8 10 12]
        @test aggregate(array1, End(), scale) == [14 16 18]
    end
    @testset "mixed locus" begin
        @test aggregate(array1, (End(), Start()), 3) == [13 16]
        @test aggregate(array1, (End(), Start()), (3, 2)) == [13 15 17]
    end
end

@testset "Aggregate with a function" begin
    @test aggregate(array1, sum, 3) == [72 99]
    @test aggregate(array1, median, 3) == [8 11]
    @test aggregate(array1, sum, (3, 2)) == [45 57 69]
end

