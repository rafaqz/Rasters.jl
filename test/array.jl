using GeoData, Test, Dates
using GeoData.LookupArrays, GeoData.Dimensions

data1 = cumsum(cumsum(ones(10, 11); dims=1); dims=2)
data2 = 2cumsum(cumsum(ones(10, 11, 1); dims=1); dims=2)
dims1 = X(10:10:100), Y(-50:10:50) 
dims2 = (dims1..., Ti([DateTime(2019)]))
refdimz = ()
mval = -9999.0
meta = NoMetadata()
nme = :test

ga2 = GeoArray(data2, dims2)
ga1 = GeoArray(data1; dims=dims1, refdims=refdimz, name=nme, metadata=meta, missingval=mval)

@test ga1 == data1
@test ga2 == data2

@testset "arary dims have been formatted" begin
    @test index(ga2) == (10.0:10:100, -50.0:10:50.0, [DateTime(2019)])
    @test dims(ga1)[1:2] == dims(ga2)[1:2]
    @test name(ga1) == :test
    @test missingval(ga1) == -9999.0
end

@testset "show" begin
    sh = sprint(show, MIME("text/plain"), ga1)
    # Test but don't lock this down too much
    @test occursin("GeoArray", sh)
    @test occursin("Y", sh)
    @test occursin("X", sh)
end
