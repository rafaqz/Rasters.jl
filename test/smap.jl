using HDF5, GeoData, Test, Statistics, Dates
using GeoData: Time, window, name
include("test_utils.jl")

# TODO example files without a login requirement
path1 = "SMAP_L4_SM_gph_20160101T223000_Vv4011_001.h5"
path2 = "SMAP_L4_SM_gph_20160102T223000_Vv4011_001.h5"

if isfile(path1)
    @testset "stack" begin
        stack = SMAPstack(path1)
        dims(stack, :soil_temp_layer1)
        keys(stack)
        names(stack)

        @testset "conversion to GeoArray" begin
            smaparray = stack["soil_temp_layer1"][Lon(1:100), Lat(1:100)]
            @test typeof(smaparray) <: GeoArray{Float32,2}
            @test size(smaparray) == (100, 100)
            @test typeof(dims(smaparray)) <: Tuple{<:Lon{<:Array{Float32,1}}, <:Lat{<:Array{Float32,1}}}
            @test typeof(refdims(smaparray)) <: Tuple{<:Time} 
            @test missingval(smaparray) == -9999.0
            @test smaparray[1] == -9999.0
            @test name(smaparray) == "soil_temp_layer1"
            # Why is tagged time different to the filename time? is that just rounded?
            @test refdims(stack) == (Time(DateTime(2016, 1, 1, 22, 28, 55, 816)),)
            @test_broken metadata(smaparray) = "not implemented yet"
        end

        @testset "conversion to GeoStack" begin
            # Stack Constructors
            # This uses too much ram! There is a lingering memory leak in HDF5.
            # stack = GeoStack(stack) 
            # keys(stack)
            # @test Symbol.(Tuple(keys(stack))) == keys(stack)
            geostack = GeoStack(stack; keys=(:baseflow_flux, :snow_mass, :soil_temp_layer1))
            keys(geostack) == (:baseflow_flux, :snow_mass, :soil_temp_layer1)
        end

        if VERSION > v"1.1-"
            @testset "copy" begin
                geoarray = zero(stack[:soil_temp_layer1])
                @test typeof(geoarray) <: GeoArray
                @test geoarray != stack[:soil_temp_layer1]
                copy!(geoarray, stack, :soil_temp_layer1)
                @test geoarray == stack[:soil_temp_layer1]
            end
        end

        @testset "window" begin
            windowedstack = SMAPstack(path1; window=(Lat(1:5), Lon(1:5), Time(1)))
            @test window(windowedstack) == (Lat(1:5), Lon(1:5), Time(1))
            windowedarray = windowedstack[:soil_temp_layer1]
            @test size(windowedarray) == (5, 5)
            @test windowedarray[1:3, 2:2] == reshape([-9999.0, -9999.0, -9999.0], 3, 1)
            @test windowedarray[1:3, 2] == [-9999.0, -9999.0, -9999.0]
            @test windowedarray[1, 2] == -9999.0 
            windowedstack = SMAPstack(path1; window=(Lat(1:5), Lon(1:5), Time(1:1)))
            windowedarray = windowedstack[:soil_temp_layer1]
            @test windowedarray[1:3, 2:2, 1] == reshape([-9999.0, -9999.0, -9999.0], 3, 1)
            @test windowedarray[1:3, 2, 1] == [-9999.0, -9999.0, -9999.0]
            @test windowedarray[1, 2, 1] == -9999.0
            windowedstack = SMAPstack(path1; window=(Time(1),))
            windowedarray = windowedstack[:soil_temp_layer1]
            @test windowedarray[1:3, 2:2] == reshape([-9999.0, -9999.0, -9999.0], 3, 1)
            @test windowedarray[1:3, 2] == [-9999.0, -9999.0, -9999.0]
            @test windowedarray[1, 2] == -9999.0
        end

    end

    @testset "series" begin
        series = SMAPseries([path1, path2])
        @test typeof(series[1]) <: SMAPstack
        # Bounds are nothing - they need to be input somehow
        @test_broken first(bounds(series, Time)) == DateTime(2016, 1, 1, 22, 28, 55, 816)
    end
    
end
