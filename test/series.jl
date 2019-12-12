using GeoData, Test, Dates
using GeoData: Time, formatdims, dims

# GeoSeries from GeoArray/GeoStack components

data1 = [1 2 3 4
         5 6 7 8]
data2 = 2 * data1
data3 = 3 * data1
data4 = 4 * data1
dimz = Lon([30, 40]), Lat((-10, 20))
ga1 = GeoArray(data1, dimz)
ga2 = GeoArray(data2, dimz)
ga1a = GeoArray(data3, dimz)
ga2a = GeoArray(data4, dimz)
stack1 = GeoStack(ga1, ga2; keys=(:ga1, :ga2))
stack2 = GeoStack(ga1a, ga2a; keys=(:ga1, :ga2))
dates =[DateTime(2017), DateTime(2018)]
series = GeoSeries([stack1, stack2], (Time(dates),));
@test issorted(dates)

@testset "getindex returns the currect types" begin
    @test typeof(series[Time(1)]) <: GeoStack{<:NamedTuple}
    @test typeof(series[Time(1)][:ga2]) <: GeoArray{Int,2}
    @test typeof(series[Time(1)][:ga2, 1, 1]) <: Int
    @test typeof(series[Time(1)][:ga2][1, 1]) <: Int
end

@testset "getindex returns the currect results" begin
    @test series[Time<|Near<|DateTime(2017)][:ga1][Lon(1), Lat(3)] === 3
    @test series[Time<|At<|DateTime(2017)][:ga1, Lon<|1, Lat<|3] === 3
    @test series[Time<|At<|DateTime(2018)][:ga2][Lon(2), Lat(4)] === 32
    @test series[Time<|At<|DateTime(2018)][:ga2, Lon<|2, Lat<|4] === 32
    @test series[Time(1)][:ga1, Lon(1), Lat(2)] == 2
    @test series[Time(1)][:ga2, Lon(2), Lat(3:4)] == [14, 16] 
end

@testset "getindex is type stable all the way down" begin
    @inferred series[Time<|At(DateTime(2017))][:ga1, Lon(1), Lat(2)]
    @inferred series[Time(1)][:ga1][Lon(1), Lat(2)]
    @inferred series[Time(1)][:ga1, Lon(1), Lat(2:4)]
    @inferred series[Time(1)][:ga1][Lon(1), Lat(2:4)]
    @inferred series[1][:ga1, Lon(1:2), Lat(:)]
    @inferred series[1][:ga1][Lon(1:2), Lat(:)]
end


# @testset "lazy view windows" begin
#     dimz = (Time<|[DateTime(2017), DateTime(2018)],)
#     dat = [stack1, stack2]
#     windowdimz = Lon(1:2), Lat(3:4)
#     series = GeoSeries(dat, dimz; window=windowdimz)
#     @test window(series) == windowdimz
#     stack = series[1]
#     @test window(stack) == windowdimz
#     @test stack[:ga1] == [3 4; 7 8]
#     @test stack[:ga1, 1, 2] == 4
# end
