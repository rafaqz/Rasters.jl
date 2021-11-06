using Rasters, Test, Dates
using Rasters.LookupArrays, Rasters.Dimensions

data1 = cumsum(cumsum(ones(10, 11); dims=1); dims=2)
data2 = 2cumsum(cumsum(ones(10, 11, 1); dims=1); dims=2)
dims1 = X(10:10:100), Y(-50:10:50) 
dims2 = (dims1..., Ti([DateTime(2019)]))
refdimz = ()
mval = -9999.0
meta = NoMetadata()
nme = :test

ga2 = Raster(data2, dims2)
ga1 = Raster(data1; dims=dims1, refdims=refdimz, name=nme, metadata=meta, missingval=mval)

@test ga1 == data1
@test ga2 == data2

@testset "arary dims have been formatted" begin
    @test index(ga2) == (10.0:10:100, -50.0:10:50.0, [DateTime(2019)])
    @test dims(ga1)[1:2] == dims(ga2)[1:2]
    @test name(ga1) == :test
    @test missingval(ga1) == -9999.0
end

@testset "constructors" begin
    x = X(Projected(1.0:1.0:2.0; crs=ProjString("+proj=")))
    y = Y(Projected(1.0:1.0:2.0; crs=ProjString("+proj=")))

    A1 = rand(Float32, x, y)
    @test A1 isa Raster{Float32,2}
    @test size(A1) == (2, 2)

    A2 = zeros(Int, x, y)
    @test A2 isa Raster{Int,2}
    @test A2 == [0 0; 0 0] 

    A3 = ones(x, y)
    @test A3 isa Raster{Float64,2}
    @test A3 == [1.0 1.0; 1.0 1.0] 

    A4 = fill(99.0, x, y)
    @test A4 isa Raster{Float64,2}
    @test A4 == [99.0 99.0; 99.0 99.0] 

    A5 = trues(x, y)
    @test A5 isa Raster{Bool,2}
    @test A5 == [true true; true true] 

    A6 = falses(x, y)
    @test A6 isa Raster{Bool,2}
    @test A6 == [false false; false false] 
end

@testset "show" begin
    sh = sprint(show, MIME("text/plain"), ga1)
    # Test but don't lock this down too much
    @test occursin("Raster", sh)
    @test occursin("Y", sh)
    @test occursin("X", sh)
end
