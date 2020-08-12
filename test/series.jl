using GeoData, Test, Dates
using GeoData: formatdims, childtype

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
series = GeoSeries([stack1, stack2], (Ti(dates),))
@test issorted(dates)

@testset "getindex returns the currect types" begin
    @test series[Ti(1)] isa GeoStack{<:NamedTuple}
    @test series[Ti(1)][:ga2] isa GeoArray{Int,2}
    @test series[Ti(1)][:ga2, 1, 1] isa Int
    @test series[Ti(1)][:ga2][1, 1] isa Int
end

@testset "properties" begin
    @test childtype(series) === GeoStack
    @test refdims(series) === ()
    # Should these be real fields? what is the use-case?  @test metadata(series) === nothing @test name(series) === ""
    @test label(series) === ""
end

@testset "getindex returns the currect results" begin
    @test series[Ti(Near(DateTime(2017)))][:ga1][Lon(1), Lat(3)] === 3
    @test series[Ti(At(DateTime(2017)))][:ga1, Lon(1), Lat(3)] === 3
    @test series[Ti(At(DateTime(2018)))][:ga2][Lon(2), Lat(4)] === 32
    @test series[Ti(At(DateTime(2018)))][:ga2, Lon(2), Lat(4)] === 32
    @test series[Ti(1)][:ga1, Lon(1), Lat(2)] == 2
    @test series[Ti(1)][:ga2, Lon(2), Lat(3:4)] == [14, 16] 
end

@testset "getindex is type stable all the way down" begin
    # @inferred series[Ti<|At(DateTime(2017))][:ga1, Lon(1), Lat(2)]
    @inferred series[Ti(1)][:ga1][Lon(1), Lat(2)]
    # @inferred series[Ti(1)][:ga1, Lon(1), Lat(2:4)]
    @inferred series[Ti(1)][:ga1][Lon(1), Lat(2:4)]
    # @inferred series[1][:ga1, Lon(1:2), Lat(:)]
    @inferred series[1][:ga1][Lon(1:2), Lat(:)]
end

@testset "lazy view windows" begin
    dimz = (Ti<|[DateTime(2017), DateTime(2018)],)
    dat = [stack1, stack2]
    window_ = Lon(1:2), Lat(3:4)
    series = GeoSeries(dat, dimz; childkwargs=(window=window_,))
    stack = series[1]
    @test GeoData.window(stack) == window_
    @test stack[:ga1] == [3 4; 7 8]
    @test stack[:ga1, 1, 2] == 4
end

@testset "setindex!" begin
    series[1] = series[1]
    @test typeof(series[1]) == typeof(series[2]) == eltype(parent(series))
    typeof(parent(series)[1]) == eltype(series)
    series[Ti(1)] = series[Ti(2)]
    @test series[Ti(1)] == series[Ti(2)]
end
