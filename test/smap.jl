# TODO example files without a login requirement
smapfile1 = "SMAP_L4_SM_gph_20160101T223000_Vv4011_001.h5"
smapfile2 = "SMAP_L4_SM_gph_20160102T223000_Vv4011_001.h5"

@testset "smapstack" begin
    smapstack = SMAPstack(smapfile1);
    dims(smapstack, "soil_temp_layer1")
    keys(smapstack)
    names(smapstack)

    @testset "conversion to GeoArray" begin
        geoarray = smapstack["soil_temp_layer1", Lon(1:100), Lat(1:100)]
        @test typeof(geoarray) <: GeoArray{Float32,2}
        @test size(geoarray) == (100, 100)
        @test typeof(dims(geoarray)) <: Tuple{<:Lon{<:Array{Float32,1}}, <:Lat{<:Array{Float32,1}}}
    end

    @testset "conversion to GeoStack" begin
        # Stack Constructors
        # This uses too much ram! There is a lingering memory leak in HDF5.
        # stack = GeoStack(smapstack) 
        # keys(stack)
        # @test Symbol.(Tuple(keys(smapstack))) == keys(smapstack)
        stack = GeoStack(smapstack; keys=(:baseflow_flux, :snow_mass, :soil_temp_layer1))
        keys(stack) == (:baseflow_flux, :snow_mass, :soil_temp_layer1)
        # stack = GeoStack(smapstack; keys=("baseflow_flux", "snow_mass", "soil_temp_layer1"))
        # keys(stack) == (:baseflow_flux, :snow_mass, :soil_temp_layer1)
    end

    @testset "copy" begin
        array = GeoArray(smapstack[:soil_temp_layer1])
        @test typeof(array) <: GeoArray
        @test array == smapstack[:soil_temp_layer1]
        copy!(array, smapstack, :soil_temp_layer2)
        @test array == smapstack[:soil_temp_layer2]
    end
end
