using GeoData, Test, Statistics, Dates
import ArchGDAL, NCDatasets, HDF5
using GeoData: Time, window, name

testpath = joinpath(dirname(pathof(GeoData)), "../test/")
include(joinpath(testpath, "test_utils.jl"))

# TODO example files without a login requirement
path1 = joinpath(testpath, "data/SMAP_L4_SM_gph_20160101T223000_Vv4011_001.h5")
path2 = joinpath(testpath, "data/SMAP_L4_SM_gph_20160102T223000_Vv4011_001.h5")

if isfile(path1) && isfile(path2)
    @testset "stack" begin
        smapstack = stack(path1)

        @testset "read" begin
            st = read(smapstack)
            @test st isa GeoStack
            @test st.data isa NamedTuple
            @test first(st.data) isa GeoArray
            @test parent(first(st.data)) isa Array
        end

        @testset "conversion to GeoArray" begin
            smaparray = smapstack["soil_temp_layer1"][Y(), X()]
            @test smaparray isa GeoArray{Float32,2}
            @test dims(smaparray) isa Tuple{<:X{<:Array{Float32,1}}, <:Y{<:Array{Float32,1}}}
            @test span(smaparray) isa Tuple{Irregular{Tuple{Float32,Float32}},Irregular{Tuple{Float32,Float32}}}
            @test span(smaparray) == (Irregular((-180.0f0, 180.0f0)), Irregular((-85.04456f0, 85.04456f0)))
            @test bounds(smaparray) == ((-180.0f0, 180.0f0), (-85.04456f0, 85.04456f0))
            @test index(smaparray) isa Tuple{Vector{Float32},Vector{Float32}}
            @test refdims(smaparray) isa Tuple{<:Ti}
            @test missingval(smaparray) == -9999.0
            @test smaparray[1] == -9999.0
            @test name(smaparray) == :soil_temp_layer1
            dt = DateTime(2016, 1, 1, 22, 30)
            step_ = Hour(3)
            @test refdims(smapstack) ==
                (Ti(dt:step_:dt; mode=Sampled(Ordered(), Regular(step_), Intervals(Start()))),)
            # Currently empty
            @test metadata(smaparray) isa Metadata{:SMAP}
        end

        @testset "conversion to GeoStack" begin
            # Stack Constructors
            # This uses too much ram! There is a lingering memory leak in HDF5.
            # geostack = GeoStack(stack)
            # @test Symbol.(Tuple(keys(stack))) == keys(geostack)
            geostack = GeoStack(smapstack; keys=(:baseflow_flux, :snow_mass, :soil_temp_layer1))
            @test keys(geostack) == (:baseflow_flux, :snow_mass, :soil_temp_layer1)
        end

        if VERSION > v"1.1-"
            @testset "copy" begin
                geoA = zero(smapstack[:soil_temp_layer1])
                @test geoA isa GeoArray
                @test geoA != smapstack[:soil_temp_layer1]
                copy!(geoA, smapstack, :soil_temp_layer1)
                @test geoA == smapstack[:soil_temp_layer1]
            end
        end

        @testset "window" begin
            windowedstack = SMAPstack(path1; window=(Y(1:5), X(1:5), Ti(1)))
            @test window(windowedstack) == (Y(1:5), X(1:5), Ti(1))
            windowedarray = windowedstack[:soil_temp_layer1]
            @test size(windowedarray) == (5, 5)
            @test windowedarray[1:3, 2:2] == reshape([-9999.0, -9999.0, -9999.0], 3, 1)
            @test windowedarray[1:3, 2] == [-9999.0, -9999.0, -9999.0]
            @test windowedarray[1, 2] == -9999.0
            windowedstack = SMAPstack(path1; window=(Y(1:5), X(1:5), Ti(1:1)))
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

        @testset "show" begin
            sh1 = sprint(show, smapstack[:soil_temp_layer1])
            # Test but don't lock this down too much
            @test occursin("GeoArray", sh1)
            @test occursin("Y", sh1)
            @test occursin("X", sh1)
            @test occursin("Time", sh1)
            sh2 = sprint(show, smapstack[:soil_temp_layer1][Y(Between(0, 100)), X(Between(1, 100))])
            # Test but don't lock this down too much
            @test occursin("GeoArray", sh2)
            @test occursin("Y", sh2)
            @test occursin("X", sh2)
            @test occursin("Time", sh2)
        end

    end

    @testset "series" begin
        smapseries = SMAPseries([path1, path2]);
        val.(dims(smapseries))
        @test smapseries[1] isa SMAPstack
        @test first(bounds(smapseries, Ti)) == DateTime(2016, 1, 1, 22, 30)
        @test last(bounds(smapseries, Ti)) == DateTime(2016, 1, 3, 1, 30)
        modified_series = modify(Array, smapseries)
        stackkeys = keys(modified_series[1])
        @test typeof(modified_series) <: GeoSeries{<:GeoStack{<:NamedTuple{stackkeys,<:Tuple{<:GeoArray{Float32,2,<:Tuple,<:Tuple,<:Array{Float32,2}},Vararg}}}}

        @testset "read" begin
            geoseries = read(smapseries)
            @test geoseries isa GeoSeries{<:GeoStack}
            @test geoseries.data isa Vector{<:GeoStack}
            @test first(geoseries.data[1].data) isa GeoArray 
        end
    end
end
