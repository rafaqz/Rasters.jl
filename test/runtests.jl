using GeoData, Test, Statistics, Dates
using NCDatasets, HDF5, ArchGDAL
using GeoData: Time, formatdims, data

data1 = cumsum(cumsum(ones(10, 11); dims=1); dims=2)
data2 = 2cumsum(cumsum(ones(10, 11, 1); dims=1); dims=2)
dims1 = Lon<|(10, 100), Lat<|(-50, 50) 
dims2 = (dims1..., Time<|[DateTime(2019)])
refdimz = ()
mval = -9999.0
meta = nothing

ga1 = GeoArray(data1, dims1, refdimz, meta, mval)
ga2 = GeoArray(data2, dims2)
# Dims have been formatted
@test val.(dims(ga2)) == val.((Lon<|LinRange(10.0, 100.0, 10), 
                               Lat<|LinRange(-50.0, 50.0, 11),
                               Time<|[DateTime(2019)]))
@test dims(ga1)[1:2] == dims(ga2)[1:2]


# GeoStack. This also tests GeoArray for the same methods.

# Constructor
stack = GeoStack(ga1, ga2; keys=(:ga1, :ga2), dims=dims)
@test typeof(parent(stack)) <: NamedTuple
@test length(parent(stack)) == 2
@test stack[:ga1] === ga1
@test stack[:ga2] === ga2

# getters
# @test data(stack, :ga1) == ga1
@test dims(stack) == formatdims(data1, dims1)
@test dims(stack, :ga1) == formatdims(data1, dims1)
@test refdims(stack) == ()
@test metadata(stack) == nothing
@test metadata(stack, :ga1) == nothing
@test Base.keys(stack) == (:ga1, :ga2)
@test values(stack) == (ga1, ga2)

# Indexing the stack is the same as indexing its child array
a = stack[:ga1][Lon<|2:4, Lat<|5:6]
@test a == stack[:ga1, Lon<|2:4, Lat<|5:6]

@inferred stack[:ga1][Lon<|2:4, Lat<|5:6] 
# FIXME: This isn't inferred, the constants don't propagate like they 
# do in the above call. Probably due to the anonymous wrapper funcion. 
# It's not actually needed for GeoStack so could be worked around
# But its also means it wont infer for other types, so should be fixed at
# the source.
@test_broken @inferred stack[:ga1, Lon<|2:4, Lat<|5:6] 

# Getindex for a whole stack of new GeoArrays
a = stack[Lon<|2:4, Lat<|5:6]
@test typeof(a) <: GeoStack
@test typeof(a[:ga1]) <: GeoArray
@test typeof(parent(a[:ga1])) <: Array
@test a[:ga1] == data1[2:4, 5:6]
@test a[:ga2] == data2[2:4, 5:6, 1:1]

# Select new arrays for the whole stack
s = select(stack, Lat<|At(-10, 10.0), Time<|DateTime(2019))
select(stack, Lat<|At(-10, 10.0), Time<|DateTime(2019))
typeof(s)
keys(s)
@test typeof(s) <: GeoStack
@test typeof(s[:ga1]) <: GeoArray
@test typeof(parent(s[:ga1])) <: Array
@test s[:ga1] == data1[:, 5:7]
@test s[:ga2] == data2[:, 5:7, 1]
@test dims(s[:ga2]) == (Lon<|LinRange(10.0, 100.0, 10), 
                        Lat<|LinRange(-10.0, 10.0, 3))
@test dims(s, :ga2) == dims(s[:ga2])
@test refdims(s[:ga2]) == (Time<|DateTime(2019),)
# @test missingval(s, :ga2) == missingval(s[:ga2])

# Select views of arrays for the whole stack
sv = selectview(stack, Lat<|Between(-4.0, 27.0), Time<|DateTime(2019))
@test typeof(sv) <: GeoStack
@test typeof(sv[:ga1]) <: GeoArray
@test typeof(parent(sv[:ga1])) <: SubArray
@test sv[:ga1] == data1[:, 6:8]
@test sv[:ga2] == data2[:, 6:8, 1]
@test dims(sv[:ga2]) == (Lon<|LinRange(10.0, 100.0, 10), Lat<|LinRange(0.0, 20.0, 3))
@test refdims(sv[:ga2]) == (Time<|DateTime(2019),)

# Stack of view-based GeoArrays
v = view(stack, Lon(2:4), Lat(5:6))
# @inferred view(stack, Lon(2:4), Lat(5:6))
@test typeof(v) <: GeoStack
@test typeof(v[:ga1]) <: GeoArray
@test typeof(parent(v[:ga1])) <: SubArray
@test v[:ga1] == view(data1, 2:4, 5:6)
@test v[:ga2] == view(data2, 2:4, 5:6, 1:1)

# New stack with specific key(s)
x = GeoStack(stack, :ga2)
@test keys(x) == (:ga2,)
y = GeoStack(stack, (:ga1, :ga2))
@test keys(y) == (:ga1, :ga2)


# GeoSeries

data1 = [1 2 3 4
         5 6 7 8]
data2 = 2data1
data3 = 3data1
data4 = 4data1
dimz = Lon<|[30, 40], Lat<|(-10, 20)
ga1 = GeoArray(data1, dimz)
ga2 = GeoArray(data2, dimz)
ga1a = GeoArray(data3, dimz)
ga2a = GeoArray(data4, dimz)
stack1 = GeoStack(ga1, ga2; keys=(:ga1, :ga2), dims=dims)
stack2 = GeoStack(ga1a, ga2a; keys=(:ga1, :ga2), dims=dims)
series = GeoSeries([stack1, stack2], (Time<|[DateTime(2017), DateTime(2018)],));

@test select(series, Time<|Near<|DateTime(2017))[:ga1][Lon(1), Lat(3)] === 3
@test select(series, Time<|At<|DateTime(2018))[:ga2][Lon(2), Lat(4)] === 32

@test select(series, Time<|DateTime(2017))[:ga1, Lon<|1, Lat<|3] === 3
@test select(series, Time<|DateTime(2018))[:ga2, Lon<|2, Lat<|4] === 32

@test series[Time(1)][:ga1, Lon(1), Lat(2)] == 2
@test series[Time(1)][:ga2, Lon(2), Lat(3:4)] == [14, 16] 
# The contained arrays have the same dimensions, this should be type stable
@inferred series[Time(1)][:ga1, Lon(1), Lat(2)]
@inferred series[Time(1)][:ga1, Lon(1), Lat(2:4)]
@inferred series[1][:ga1, Lon(1:2), Lat(:)]


# External sources

geturl(url) = begin
    fname = splitdir(url)[2]
    isfile(fname) || download(url, fname)
    fname
end

# NCDatasets

ncexamples = "https://www.unidata.ucar.edu/software/netcdf/examples/"
ncsingle = joinpath(ncexamples, "tos_O1_2001-2002.nc")
ncmulti = joinpath(ncexamples, "test_echam_spectral.nc")

ncarray = NCarray(geturl(ncsingle))
@test typeof(ncarray) <: NCarray{Union{Missing,Float32},3}
@test typeof(ncarray[Time(1)]) <: GeoArray{Union{Missing,Float32},2}
@test typeof(dims(ncarray)) <: Tuple{<:Lon,<:Lat,<:Time}
@test length(val(dims(dims(ncarray), Time))) == 24

ncstack = NCstack(geturl(ncmulti))
@test typeof(ncstack) <: NCstack{String}
# Rebuilds as a regular GeoArray
@test typeof(ncstack[:albedo]) <: NCarray{Union{Missing,Float32},3} 
@test typeof(ncstack[:albedo, 2, 3, 1]) <: Float32
@test typeof(ncstack[:albedo, :, 3, 1]) <: GeoArray{Union{Missing,Float32},1}
@test typeof(dims(ncstack, :albedo)) <: Tuple{<:Lon,<:Lat,<:Time}
@test length(keys(ncstack)) == 136
@test first(keys(ncstack)) == "abso4"
# Test some DimensionalData.jl tools work
# Time dim should be reduced to length 1 by mean
@test axes(mean(ncstack[:albedo, Lat(1:20)], dims=Time)) == (Base.OneTo(192), Base.OneTo(20), Base.OneTo(1))
ncstack[:albedo, Lat(1:20)]
metadata(ncstack, :albedo)
keys(ncstack)
array = ncstack[:albedo][Time(4:6), Lon(1), Lat(2)] 
@test array == ncstack[:albedo, Time(4:6), Lon(1), Lat(2)] 
size(array) == (3,)

ncmultistack = NCstack([geturl(ncsingle)])
ncmultistack = NCstack((geturl(ncsingle),))
@test typeof(ncmultistack[:tos]) <: NCarray{Union{Missing,Float32},3}
@test typeof(ncstack[:albedo, 2, 3, 1]) <: Float32
@test typeof(ncstack[:albedo, Time(1)]) <: GeoArray{Union{Missing,Float32},2}

geoarray = GeoArray(ncarray)
@test typeof(dims(geoarray)) <: Tuple{<:Lon,<:Lat,<:Time}
@test dims(geoarray) == dims(ncarray)
@test metadata(geoarray) == metadata(ncarray)
@test ismissing(missingval(geoarray))


# ArchGDAL

gdal_url = "https://download.osgeo.org/geotiff/samples/gdal_eg/cea.tif"
gdalarray = GDALarray(geturl(gdal_url))
@test ndims(gdalarray) == 3
@test size(gdalarray) == (514, 515, 1)
@test typeof(gdalarray) <: GDALarray{UInt8,3}
@test typeof(dims(gdalarray)) <: Tuple{Lon,<:Lat,<:Band}
@test length(val(dims(dims(gdalarray), Lon))) == 514

# getindex to memory-backed GeoArray
@time gdalarray[Lon(1:50), Lat(1), Band(1)]
# Doesn't handle returning a single value
@test_broken typeof(gdalarray[Lon(1), Lat(1), Band(1)]) <: UInt8
@test_broken typeof(gdalarray[1, 1, 1]) <: UInt8
@time geoarray = gdalarray[Lon(1), Lat(1:50), Band(1)]

@time typeof(geoarray) <: GeoArray{UInt8,1} 
@test typeof(dims(geoarray)) <: Tuple{<:Lat}
@test typeof(refdims(geoarray)) <: Tuple{<:Lon,<:Band} 
@test metadata(geoarray) == metadata(gdalarray)
@test missingval(geoarray) == -1.0e10

gdalstack = GDALstack((a=geturl(gdal_url),))
parent(gdalstack)

# Broken: gdal read() indexing is non-standard
# band is first when it is last in the returned array.
# and at least one unitrange is required for dispatch.
@test_broken gdalstack[:a][Lat([2,3]), Lon(1), Band(1)]
@test_broken gdalstack[:a][Lat(1), Lon(1), Band(1)]


# SMAP

smapfile = "SMAP_L4_SM_gph_20160101T223000_Vv4011_001.h5"
smapstack = SMAPstack(smapfile);
dims(smapstack, "soil_temp_layer1")
keys(smapstack)
names(smapstack)
geoarray = smapstack["soil_temp_layer1", Lon(1:100), Lat(1:100)]
geoarray == smapstack["soil_temp_layer1"][Lon(1:100), Lat(1:100)]
@test typeof(geoarray) <: GeoArray{Float32,2}
@test size(geoarray) == (100, 100)
@test typeof(dims(geoarray)) <: Tuple{<:Lon{<:Array{Float32,1}}, <:Lat{<:Array{Float32,1}}}


