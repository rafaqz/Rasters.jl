using NCDatasets, GeoData, Test, Statistics, Dates, CFTime
using GeoData: Time, window, name
include("test_utils.jl")

ncexamples = "https://www.unidata.ucar.edu/software/netcdf/examples/"
ncsingle = geturl(joinpath(ncexamples, "tos_O1_2001-2002.nc"))
ncmulti = geturl(joinpath(ncexamples, "test_echam_spectral.nc"))

@testset "NCarray" begin
    ncarray = NCDarray(ncsingle)

    @testset "array properties" begin
        @test size(ncarray) == (180, 170, 24)
        @test typeof(ncarray) <: NCDarray
    end

    @testset "dimensions" begin
        @test ndims(ncarray) == 3
        @test length.(val.(dims(ncarray))) == (180, 170, 24)
        @test typeof(dims(ncarray)) <: Tuple{<:Lon,<:Lat,<:Time}
        @test refdims(ncarray) == ()
        @test bounds(ncarray) == ((1.0, 359.0), (-79.5, 89.5), (DateTime360Day(2001, 1, 16), DateTime360Day(2002, 12, 16)))
    end

    @testset "other fields" begin
        @test window(ncarray) == ()
        @test ismissing(missingval(ncarray))
        @test typeof(metadata(ncarray)) <: NCDmetadata # TODO make this a namedtuple
        @test name(ncarray) == "tos"
    end

    @testset "indexing" begin
        @test typeof(ncarray[Time(1)]) <: GeoArray{Union{Missing,Float32},2}
        @test typeof(ncarray[Lat(1), Time(1)]) <: GeoArray{Union{Missing,Float32},1}
        @test typeof(ncarray[Lon(1), Time(1)]) <: GeoArray{Union{Missing,Float32},1}
        @test typeof(ncarray[Lon(1), Lat(1), Time(1)]) <: Missing
        @test typeof(ncarray[Lon(30), Lat(30), Time(1)]) <: Float32
        @test typeof(ncarray[30, 30, 2]) <: Float32
    end

    @testset "selectors" begin
        a = ncarray[Lon(At(21.0)), Lat(Between(50, 52)), Time(Near(DateTime360Day(2002, 12)))]
        @test bounds(a) == ((50.5, 51.5),)
        x = ncarray[Lon(Near(150)), Lat(Near(30)), Time(1)]
        @test typeof(x) <: Float32
        # TODO make sure we are getting the right cell.
    end

    @testset "conversion to GeoArray" begin
        geoarray = ncarray[Lon(1:50), Lat(20:20), Time(1)]
        @test size(geoarray) == (50, 1)
        @test eltype(geoarray) <: Union{Missing,Float32}
        @time typeof(geoarray) <: GeoArray{Float32,1} 
        @test typeof(dims(geoarray)) <: Tuple{<:Lon,<:Lat}
        @test typeof(refdims(geoarray)) <: Tuple{<:Time} 
        @test metadata(geoarray) == metadata(ncarray)
        @test ismissing(missingval(geoarray))
        @test name(geoarray) == "tos"
    end

    @testset "window" begin
        ds = Dataset(ncsingle)
        ds["tos"][101:105, 51:55, 1][1:3, 2:2]
        windowedarray = NCDarray(ncsingle; window=(Lat(51:55), Lon(101:105), Time(1)))
        @test size(windowedarray) == (5, 5)
        @test window(windowedarray) == (101:105, 51:55, 1)
        @test ndims(windowedarray) == 2
        @test windowedarray[1:3, 2:2] == reshape([297.3289f0, 297.44012f0, 297.4756f0], 3, 1)
        @test windowedarray[1:3, 2] == [297.3289f0, 297.44012f0, 297.4756f0]
        @test windowedarray[1, 2] == 297.3289f0
    end

    @testset "save" begin
        geoarray = GeoArray(ncarray)
        @test size(geoarray) == size(ncarray)
        filename = tempname()
        GeoData.write(filename, NCDarray, geoarray)
        saved = GeoArray(NCDarray(filename))
        @test size(saved) == size(geoarray)
        @test refdims(saved) == refdims(geoarray)
        @test missingval(saved) === missingval(geoarray)
        @test_broken metadata(saved) == metadata(geoarray)
        @test GeoData.name(saved) == GeoData.name(geoarray)
        @test_broken all(metadata.(dims(saved)) .== metadata.(dims(geoarray)))
        @test all(DimensionalData.grid.(dims(saved)) .== DimensionalData.grid.(dims(geoarray)))
        @test_broken typeof(dims(saved)) == typeof(dims(geoarray))
        @test val(dims(saved)[3]) == val(dims(geoarray)[3])
        @test all(val.(dims(saved)) .== val.(dims(geoarray)))
        @test all(parent(saved) .=== parent(geoarray))
        @test_broken typeof(saved) == typeof(geoarray)
    end

end

@testset "NCDstack" begin
    ncstack = NCDstack(ncmulti)

    @testset "load ncstack" begin
        @test typeof(ncstack) <: NCDstack{String}
        @test ismissing(missingval(ncstack))
        @test typeof(metadata(ncstack)) <: NCDmetadata
        @test refdims(ncstack) == ()
        # Loads child as a regular GeoArray
        @test typeof(ncstack[:albedo]) <: GeoArray{Union{Missing,Float32},3} 
        @test typeof(ncstack[:albedo, 2, 3, 1]) <: Float32
        @test typeof(ncstack[:albedo, :, 3, 1]) <: GeoArray{Union{Missing,Float32},1}
        @test typeof(dims(ncstack, :albedo)) <: Tuple{<:Lon,<:Lat,<:Time}
        @test typeof(keys(ncstack)) == NTuple{131,Symbol}
        @test first(keys(ncstack)) == :abso4
        @test typeof(metadata(ncstack, :albedo)) <: NCDmetadata
        @test metadata(ncstack, :albedo)["institution"] == "Max-Planck-Institute for Meteorology"
        # Test some DimensionalData.jl tools work
        # Time dim should be reduced to length 1 by mean
        @test axes(mean(ncstack[:albedo, Lat(1:20)] , dims=GeoData.Time)) == 
              (Base.OneTo(192), Base.OneTo(20), Base.OneTo(1))
        geoarray = ncstack[:albedo][Time(4:6), Lon(1), Lat(2)] 
        @test geoarray == ncstack[:albedo, Time(4:6), Lon(1), Lat(2)] 
        @test size(geoarray) == (3,)
    end

    if VERSION > v"1.1-"
        @testset "copy" begin
            geoarray = GeoArray(ncstack[:albedo])
            copy!(geoarray, ncstack, :albedo)
        end
    end

    @testset "indexing" begin
        ncmultistack = NCDstack([geturl(ncsingle)])
        ncmultistack = NCDstack((geturl(ncsingle),))
        @test typeof(dims(ncmultistack)) <: Tuple{<:Lon,<:Lat,<:Time}
        @test typeof(ncmultistack[:tos]) <: GeoArray{Union{Missing,Float32},3}
        @test typeof(ncmultistack[:tos, Time(1)]) <: GeoArray{Union{Missing,Float32},2}
        @test typeof(ncmultistack[:tos, Lat(1), Time(1)]) <: GeoArray{Union{Missing,Float32},1}
        @test typeof(ncmultistack[:tos, 8, 30, 10]) <: Float32
    end

    @testset "window" begin
        windowedstack = NCDstack(ncmulti; window=(Lat(1:5), Lon(1:5), Time(1)))
        @test window(windowedstack) == (Lat(1:5), Lon(1:5), Time(1))
        windowedarray = windowedstack[:albedo]
        @test size(windowedarray) == (5, 5)
        @test windowedarray[1:3, 2:2] == reshape([0.84936917f0, 0.8776228f0, 0.87498736f0], 3, 1)
        @test windowedarray[1:3, 2] == [0.84936917f0, 0.8776228f0, 0.87498736f0]
        @test windowedarray[1, 2] == 0.84936917f0
        windowedstack = NCDstack(ncmulti; window=(Lat(1:5), Lon(1:5), Time(1:1)))
        windowedarray = windowedstack[:albedo]
        @test windowedarray[1:3, 2:2, 1:1] == reshape([0.84936917f0, 0.8776228f0, 0.87498736f0], 3, 1, 1)
        @test windowedarray[1:3, 2:2, 1] == reshape([0.84936917f0, 0.8776228f0, 0.87498736f0], 3, 1)
        @test windowedarray[1:3, 2, 1] == [0.84936917f0, 0.8776228f0, 0.87498736f0]
        @test windowedarray[1, 2, 1] == 0.84936917f0 
        windowedstack = NCDstack(ncmulti; window=(Time(1),))
        windowedarray = windowedstack[:albedo]
        @test windowedarray[1:3, 2:2] == reshape([0.84936917f0, 0.8776228f0, 0.87498736f0], 3, 1)
        @test windowedarray[1:3, 2] == [0.84936917f0, 0.8776228f0, 0.87498736f0]
        @test windowedarray[1, 2] ==  0.84936917f0
    end

    @testset "conversion to GeoStack" begin
        geostack = GeoStack(ncstack)
        @test Symbol.(Tuple(keys(geostack))) == keys(ncstack)
        smallstack = GeoStack(ncstack; keys=(:albedo, :evap, :runoff))
        @test keys(smallstack) == (:albedo, :evap, :runoff)
    end

    @testset "save" begin
        geostack = GeoStack(ncstack);
        filename = tempname()
        GeoData.write(filename, NCDstack, geostack)
        saved = GeoStack(NCDstack(filename))
        @test keys(saved) == keys(geostack)
        @test metadata(saved) == metadata(geostack)
        @test first(values(saved)) == first(values(geostack))
    end

end


