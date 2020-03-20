using NCDatasets, GeoData, Test, Statistics, Dates, CFTime
using GeoData: window, name
include("test_utils.jl")

ncexamples = "https://www.unidata.ucar.edu/software/netcdf/examples/"
ncsingle = geturl(joinpath(ncexamples, "tos_O1_2001-2002.nc"))
ncmulti = geturl(joinpath(ncexamples, "test_echam_spectral.nc"))

@testset "NCarray" begin
    ncarray = NCDarray(ncsingle)

    @testset "array properties" begin
        @test size(ncarray) == (180, 170, 24)
        @test ncarray isa NCDarray
    end

    @testset "dimensions" begin
        @test ndims(ncarray) == 3
        @test length.(val.(dims(ncarray))) == (180, 170, 24)
        @test dims(ncarray) isa Tuple{<:Lon,<:Lat,<:Ti}
        @test refdims(ncarray) == ()
        @test bounds(ncarray) == ((1.0, 359.0), (-79.5, 89.5), (DateTime360Day(2001, 1, 16), DateTime360Day(2002, 12, 16)))
    end

    @testset "other fields" begin
        @test window(ncarray) == ()
        @test ismissing(missingval(ncarray))
        @test metadata(ncarray) isa NCDarrayMetadata # TODO make this a namedtuple
        @test name(ncarray) == "tos"
    end

    @testset "indexing" begin
        @test ncarray[Ti(1)] isa GeoArray{<:Any,2}
        @test ncarray[Lat(1), Ti(1)] isa GeoArray{<:Any,1}
        @test ncarray[Lon(1), Ti(1)] isa GeoArray{<:Any,1}
        @test ncarray[Lon(1), Lat(1), Ti(1)] isa Missing
        @test ncarray[Lon(30), Lat(30), Ti(1)] isa Float32
        @test ncarray[30, 30, 2] isa Float32
    end

    @testset "selectors" begin
        a = ncarray[Lon(At(21.0)), Lat(Between(50, 52)), Ti(Near(DateTime360Day(2002, 12)))]
        @test bounds(a) == ((50.5, 51.5),)
        x = ncarray[Lon(Near(150)), Lat(Near(30)), Ti(1)]
        @test x isa Float32
        # TODO make sure we are getting the right cell.
    end

    @testset "conversion to GeoArray" begin
        geoarray = ncarray[Lon(1:50), Lat(20:20), Ti(1)]
        @test size(geoarray) == (50, 1)
        # @test eltype(geoarray) <: Union{Missing,Float32}
        @time geoarray isa GeoArray{Float32,1}
        @test dims(geoarray) isa Tuple{<:Lon,<:Lat}
        @test refdims(geoarray) isa Tuple{<:Ti}
        @test metadata(geoarray) == metadata(ncarray)
        @test ismissing(missingval(geoarray))
        @test name(geoarray) == "tos"
    end

    @testset "window" begin
        ds = Dataset(ncsingle)
        ds["tos"][101:105, 51:55, 1][1:3, 2:2]
        windowedarray = NCDarray(ncsingle; window=(Lat(51:55), Lon(101:105), Ti(1)))
        @test size(windowedarray) == (5, 5)
        @test window(windowedarray) == (101:105, 51:55, 1)
        @test ndims(windowedarray) == 2
        @test windowedarray[1:3, 2:2] == reshape([297.3289f0, 297.44012f0, 297.4756f0], 3, 1)
        @test windowedarray[1:3, 2] == [297.3289f0, 297.44012f0, 297.4756f0]
        @test windowedarray[1, 2] == 297.3289f0
    end

    @testset "save" begin
        # TODO save and load subset
        geoarray = GeoArray(ncarray)
        metadata(geoarray)
        @test size(geoarray) == size(ncarray)
        filename = tempname()
        GeoData.write(filename, NCDarray, geoarray)
        saved = GeoArray(NCDarray(filename))
        @test size(saved) == size(geoarray)
        @test refdims(saved) == refdims(geoarray)
        @test missingval(saved) === missingval(geoarray)
        @test metadata(saved) == metadata(geoarray)
        @test GeoData.name(saved) == GeoData.name(geoarray)
        @test_broken all(metadata.(dims(saved)) .== metadata.(dims(geoarray)))
        @test all(DimensionalData.grid.(dims(saved)) .== DimensionalData.grid.(dims(geoarray)))
        @test dims(saved) isa typeof(dims(geoarray))
        @test val(dims(saved)[3]) == val(dims(geoarray)[3])
        @test all(val.(dims(saved)) .== val.(dims(geoarray)))
        @test all(data(saved) .=== data(geoarray))
        @test saved isa typeof(geoarray)
    end

end

@testset "NCDstack" begin
    ncstack = NCDstack(ncmulti)

    @testset "load ncstack" begin
        @test ncstack isa NCDstack{String}
        @test ismissing(missingval(ncstack))
        @test metadata(ncstack) isa NCDstackMetadata
        @test refdims(ncstack) == ()
        # Loads child as a regular GeoArray
        @test ncstack[:albedo] isa GeoArray{<:Any,3}
        @test ncstack[:albedo, 2, 3, 1] isa Float32
        @test ncstack[:albedo, :, 3, 1] isa GeoArray{<:Any,1}
        @test dims(ncstack, :albedo) isa Tuple{<:Lon,<:Lat,<:Ti}
        @test keys(ncstack) isa NTuple{131,Symbol}
        @test first(keys(ncstack)) == :abso4
        @test metadata(ncstack, :albedo) isa NCDstackMetadata
        @test metadata(ncstack, :albedo)["institution"] == "Max-Planck-Institute for Meteorology"
        # Test some DimensionalData.jl tools work
        # Time dim should be reduced to length 1 by mean
        @test axes(mean(ncstack[:albedo, Lat(1:20)] , dims=Ti)) ==
              (Base.OneTo(192), Base.OneTo(20), Base.OneTo(1))
        geoarray = ncstack[:albedo][Ti(4:6), Lon(1), Lat(2)]
        @test geoarray == ncstack[:albedo, Ti(4:6), Lon(1), Lat(2)]
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
        @test dims(ncmultistack) isa Tuple{<:Lon,<:Lat,<:Ti}
        @test ncmultistack[:tos] isa GeoArray{<:Any,3}
        @test ncmultistack[:tos, Ti(1)] isa GeoArray{<:Any,2}
        @test ncmultistack[:tos, Lat(1), Ti(1)] isa GeoArray{<:Any,1}
        @test ncmultistack[:tos, 8, 30, 10] isa Float32
    end

    @testset "window" begin
        windowedstack = NCDstack(ncmulti; window=(Lat(1:5), Lon(1:5), Ti(1)))
        @test window(windowedstack) == (Lat(1:5), Lon(1:5), Ti(1))
        windowedarray = windowedstack[:albedo]
        @test size(windowedarray) == (5, 5)
        @test windowedarray[1:3, 2:2] == reshape([0.84936917f0, 0.8776228f0, 0.87498736f0], 3, 1)
        @test windowedarray[1:3, 2] == [0.84936917f0, 0.8776228f0, 0.87498736f0]
        @test windowedarray[1, 2] == 0.84936917f0
        windowedstack = NCDstack(ncmulti; window=(Lat(1:5), Lon(1:5), Ti(1:1)))
        windowedarray = windowedstack[:albedo]
        @test windowedarray[1:3, 2:2, 1:1] == reshape([0.84936917f0, 0.8776228f0, 0.87498736f0], 3, 1, 1)
        @test windowedarray[1:3, 2:2, 1] == reshape([0.84936917f0, 0.8776228f0, 0.87498736f0], 3, 1)
        @test windowedarray[1:3, 2, 1] == [0.84936917f0, 0.8776228f0, 0.87498736f0]
        @test windowedarray[1, 2, 1] == 0.84936917f0
        windowedstack = NCDstack(ncmulti; window=(Ti(1),))
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
        GeoData.write(filename, NCDstack, geostack);
        saved = GeoStack(NCDstack(filename))
        @test keys(saved) == keys(geostack)
        @test metadata(saved)["advection"] == "Lin & Rood"
        @test metadata(saved) == metadata(geostack)
        @test first(values(saved)) == first(values(geostack))
    end

end


