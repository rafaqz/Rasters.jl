using GeoData, Test, Statistics, Dates
using GeoData: Time, formatdims

# TODO some tests, maybe just integration and metadata.
# Everything else was moved to DimensionalData.


data1 = cumprod(cumsum(ones(10, 11); dims=1); dims=2)
data2 = cumsum(cumsum(ones(10, 11, 1); dims=1); dims=2)
dims1 = Lon((10, 100)), Lat((-50, 50)) 
dims2 = (dims1..., Time([DateTime(2019)]))
refdimz = ()
missingval = -9999.0
meta = nothing

ga1 = GeoArray(data1, dims1, refdimz, meta, missingval)
ga2 = GeoArray(data2, dims2)
# Dims have been formatted
@test val.(dims(ga2)) == val.((Lon(LinRange(10.0, 100.0, 10)), 
                               Lat(LinRange(-50.0, 50.0, 11)),
                               Time([DateTime(2019)])))
@test dims(ga1)[1:2] == dims(ga2)[1:2]


# GeoStack integration tests. This also tests GeoArray

stack = GeoStack(ga1, ga2; keys=(:ga1, :ga2), dims=dims)
@test stack[:ga1] === ga1
@test stack[:ga2] === ga2

s = select(stack, (Lat((-10, 10.0)), Time(DateTime(2019))))
@test typeof(s) <: GeoStack
@test s[:ga1] == data1[:, 5:7]
@test s[:ga2] == data2[:, 5:7, 1]
@test dims(s[:ga2]) == (Lon(LinRange(10.0, 100.0, 10)), 
                        Lat(LinRange(-10.0, 10.0, 3)))
@test refdims(s[:ga2]) == (Time(DateTime(2019)),)


v = view(stack, Lon(2:4), Lat(5:6))
@test typeof(v) <: GeoStack
@test v[:ga1] == view(data1, 2:4, 5:6)
@test v[:ga2] == view(data2, 2:4, 5:6, 1:1)

a = stack[Lon(2:4), Lat(5:6)]
@test typeof(a) <: GeoStack
@test a[:ga1] == data1[2:4, 5:6]
@test a[:ga2] == data2[2:4, 5:6, 1:1]
