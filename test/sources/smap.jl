using HDF5, GeoData, Test, Statistics, Dates, Plots
using GeoData: Time, window, name
testpath = joinpath(dirname(pathof(GeoData)), "../test/")
include(joinpath(testpath, "test_utils.jl"))

# TODO example files without a login requirement
path1 = joinpath(testpath, "data/SMAP_L4_SM_gph_20160101T223000_Vv4011_001.h5")
path2 = joinpath(testpath, "data/SMAP_L4_SM_gph_20160102T223000_Vv4011_001.h5")

@testset "stack" begin
    stack = SMAPstack(path1)

    @testset "conversion to GeoArray" begin
        smaparray = stack["soil_temp_layer1"][Lon(1:100), Lat(1:100)]
        @test smaparray isa GeoArray{Float32,2}
        @test size(smaparray) == (100, 100)
        @test dims(smaparray) isa Tuple{<:Lon{<:Array{Float32,1}}, <:Lat{<:Array{Float32,1}}}
        @test refdims(smaparray) isa Tuple{<:Ti}
        @test missingval(smaparray) == -9999.0
        @test smaparray[1] == -9999.0
        @test name(smaparray) == "soil_temp_layer1"
        # Why is tagged time different to the filename time? is that just rounded?
        dt = DateTime(2016, 1, 1, 22, 30)
        step_ = Hour(3)
        @test refdims(stack) ==
            (Ti(dt:step_:dt; mode=Sampled(Ordered(), Regular(step_), Intervals(Start()))),)
        @test_broken metadata(smaparray) = "not implemented yet"
    end

    @testset "conversion to GeoStack" begin
        # Stack Constructors
        # This uses too much ram! There is a lingering memory leak in HDF5.
        # geostack = GeoStack(stack)
        # @test Symbol.(Tuple(keys(stack))) == keys(geostack)
        geostack = GeoStack(stack; keys=(:baseflow_flux, :snow_mass, :soil_temp_layer1))
        @test keys(geostack) == (:baseflow_flux, :snow_mass, :soil_temp_layer1)
    end

    if VERSION > v"1.1-"
        @testset "copy" begin
            geoarray = zero(stack[:soil_temp_layer1])
            @test geoarray isa GeoArray
            @test geoarray != stack[:soil_temp_layer1]
            copy!(geoarray, stack, :soil_temp_layer1)
            @test geoarray == stack[:soil_temp_layer1]
        end
    end

    @testset "window" begin
        windowedstack = SMAPstack(path1; window=(Lat(1:5), Lon(1:5), Ti(1)))
        @test window(windowedstack) == (Lat(1:5), Lon(1:5), Ti(1))
        windowedarray = windowedstack[:soil_temp_layer1]
        @test size(windowedarray) == (5, 5)
        @test windowedarray[1:3, 2:2] == reshape([-9999.0, -9999.0, -9999.0], 3, 1)
        @test windowedarray[1:3, 2] == [-9999.0, -9999.0, -9999.0]
        @test windowedarray[1, 2] == -9999.0
        windowedstack = SMAPstack(path1; window=(Lat(1:5), Lon(1:5), Ti(1:1)))
        windowedarray = windowedstack[:soil_temp_layer1]
        @test windowedarray[1:3, 2:2, 1] == reshape([-9999.0, -9999.0, -9999.0], 3, 1)
        @test windowedarray[1:3, 2, 1] == [-9999.0, -9999.0, -9999.0]
        @test windowedarray[1, 2, 1] == -9999.0
        windowedstack = SMAPstack(path1; window=(Ti(1),))
        windowedarray = windowedstack[:soil_temp_layer1]
        @test windowedarray[1:3, 2:2] == reshape([-9999.0, -9999.0, -9999.0], 3, 1)
        @test windowedarray[1:3, 2] == [-9999.0, -9999.0, -9999.0]
        @test windowedarray[1, 2] == -9999.0
    end

end

@testset "series" begin
    series = SMAPseries([path1, path2]);
    val.(dims(series))
    @test series[1] isa SMAPstack
    @test first(bounds(series, Ti)) == DateTime(2016, 1, 1, 22, 30)
    @test last(bounds(series, Ti)) == DateTime(2016, 1, 3, 1, 30)
    modified_series = modify(Array, series)
    stackkeys = keys(modified_series[1])
    @test typeof(modified_series) <: GeoSeries{<:GeoStack{<:NamedTuple{stackkeys,<:Tuple{<:GeoArray{Float32,2,<:Tuple,<:Tuple,<:Array{Float32,2}},Vararg}}}}
end
