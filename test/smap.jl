# TODO example files without a login requirement
path1 = "SMAP_L4_SM_gph_20160101T223000_Vv4011_001.h5"
path2 = "SMAP_L4_SM_gph_20160102T223000_Vv4011_001.h5"

@testset "stack" begin
    stack = SMAPstack(path1);
    dims(stack, :soil_temp_layer1)
    keys(stack)
    names(stack)

    @testset "conversion to GeoArray" begin
        geoarray = stack["soil_temp_layer1", Lon(1:100), Lat(1:100)]
        dims(stack, :soil_temp_layer1)
        @test typeof(geoarray) <: GeoArray{Float32,2}
        @test size(geoarray) == (100, 100)
        @test typeof(dims(geoarray)) <: Tuple{<:Lon{<:Array{Float32,1}}, <:Lat{<:Array{Float32,1}}}
        @test typeof(refdims(geoarray)) <: Tuple{<:Time} 
        @test missingval(geoarray) == -9999.0
        @test geoarray[1] == -9999.0
        @test name(geoarray) == :soil_temp_layer1
        # Why is tagged time different to the filename time? is that just rounded?
        @test refdims(stack) == (Time(DateTime(2016, 1, 1, 22, 28, 55, 816)),)
        @test_broken metadata(geoarray) = "not implemented yet"
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

    @testset "copy" begin
        array = GeoArray(stack[:soil_temp_layer1])
        @test typeof(array) <: GeoArray
        @test array == stack[:soil_temp_layer1]
        copy!(array, stack, :soil_temp_layer2)
        @test array == stack[:soil_temp_layer2]
    end

end

@testset "series" begin
    series = SMAPseries([path1, path2])
    @test typeof(series[1]) <: SMAPstack
    @test first(bounds(series, Time)) == DateTime(2016, 1, 1, 22, 28, 55, 816)
end
