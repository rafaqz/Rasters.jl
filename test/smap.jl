# TODO example files without a login requirement
path1 = "SMAP_L4_SM_gph_20160101T223000_Vv4011_001.h5"
path2 = "SMAP_L4_SM_gph_20160102T223000_Vv4011_001.h5"

@testset "stack" begin
    stack = SMAPstack(path1);
    dims(stack, :soil_temp_layer1)
    keys(stack)
    names(stack)

    @testset "conversion to GeoArray" begin
        A = stack["soil_temp_layer1", Lon(1:100), Lat(1:100)]
        dims(stack, :soil_temp_layer1)
        @test typeof(A) <: GeoArray{Float32,2}
        @test size(A) == (100, 100)
        @test typeof(dims(A)) <: Tuple{<:Lon{<:Array{Float32,1}}, <:Lat{<:Array{Float32,1}}}
        @test typeof(refdims(A)) <: Tuple{<:Time} 
        @test missingval(A) == -9999.0
        @test A[1] == -9999.0
        @test name(A) == :soil_temp_layer1
        # Why is tagged time different to the filename time? is that just rounded?
        @test refdims(stack) == (Time(DateTime(2016, 1, 1, 22, 28, 55, 816)),)
        @test_broken metadata(A) = "not implemented yet"
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
        A = GeoArray(stack[:soil_temp_layer1])
        @test typeof(A) <: GeoArray
        @test A == stack[:soil_temp_layer1]
        copy!(A, stack, :soil_temp_layer2)
        @test A == stack[:soil_temp_layer2]
    end

end

@testset "series" begin
    series = SMAPseries([path1, path2])
    @test typeof(series[1]) <: SMAPstack
    @test first(bounds(series, Time)) == DateTime(2016, 1, 1, 22, 28, 55, 816)
end
