using GeoData, Test, Dates
using GeoData: formatdims, dims

data1 = cumsum(cumsum(ones(10, 11); dims=1); dims=2)
data2 = 2cumsum(cumsum(ones(10, 11, 1); dims=1); dims=2)
dims1 = Lon<|(10, 100), Lat<|(-50, 50) 
dims2 = (dims1..., Ti([DateTime(2019)]))
refdimz = ()
mval = -9999.0
meta = nothing
nme = :test

# Formatting only occurs in shorthand constructors
ga2 = GeoArray(data2, dims2)
ga1 = GeoArray(data1, formatdims(data1, dims1), refdimz, nme, meta, mval)

@testset "arary dims have been formatted" begin
    @test val.(dims(ga2)) == val.((Lon<|LinRange(10.0, 100.0, 10), 
                                   Lat<|LinRange(-50.0, 50.0, 11),
                                   Ti<|[DateTime(2019)]))
    @test dims(ga1)[1:2] == dims(ga2)[1:2]
end
