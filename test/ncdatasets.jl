ncexamples = "https://www.unidata.ucar.edu/software/netcdf/examples/"
ncsingle = geturl(joinpath(ncexamples, "tos_O1_2001-2002.nc"))
ncmulti = geturl(joinpath(ncexamples, "test_echam_spectral.nc"))

@testset "array" begin
    array = NCarray(ncsingle)

    @testset "array data" begin
        @test size(array) == (180, 170, 24)
        @test typeof(array) <: NCarray{Union{Missing,Float32},3}
    end

    @testset "dimensions" begin
        @test ndims(array) == 3
        @test length.(val.(dims(array))) == (180, 170, 24)
        @test typeof(dims(array)) <: Tuple{<:Lon,<:Lat,<:Time}
        @test refdims(array) == ()
        @test bounds(array) == ((1.0, 359.0), (-79.5, 89.5), (DateTime360Day(2001, 1, 16), DateTime360Day(2002, 12, 16)))
    end

    @testset "other fields" begin
        @test window(array) == ()
        @test ismissing(missingval(array))
        @test typeof(metadata(array)) <: Dict # TODO make this a namedtuple
        @test_broken metadata(array).filepath == "tos_O1_2001-2002.nc"
        @test name(array) == :tos
    end

    @testset "indexing" begin
        @test typeof(array[Time(1)]) <: GeoArray{Union{Missing,Float32},2}
        @test typeof(array[Lat(1), Band(1)]) <: GeoArray{UInt8,1}
        @test typeof(array[Lon(1), Band(1)]) <: GeoArray{UInt8,1}
        @test typeof(array[Lon(1), Lat(1), Time(1)]) <: Missing
        @test typeof(array[Lon(30), Lat(30), Time(1)]) <: Float32
        @test typeof(array[30, 30, 2]) <: Float32
    end

    @testset "selectors" begin
        a = array[Lon(At(21.0)), Lat(Between(50, 52)), Time(Near(DateTime360Day(2002, 12)))]
        @test bounds(a) == ((50.5, 51.5),)
        x = array[Lon(Near(150)), Lat(Near(30)), Time(1)]
        @test typeof(x) <: Float32
        # TODO make sure we are getting the right cell.
    end

    @testset "conversion to GeoArray" begin
        geoarray = array[Lon(1:50), Lat(20:20), Time(1)]
        @test size(geoarray) == (50, 1)
        @test eltype(geoarray) <: Union{Missing,Float32}
        @test typeof(dims(geoarray)) <: Tuple{<:Lon,<:Lat}
        @test typeof(refdims(geoarray)) <: Tuple{<:Time} 
        @test metadata(geoarray) == metadata(array)
        @test ismissing(missingval(geoarray))
        @test name(array) == :tos
    end
end

@testset "stack" begin
    stack = NCstack(geturl(ncmulti))

    @testset "load stack" begin
        @test typeof(stack) <: NCstack{String}
        @test ismissing(missingval(stack))
        @test typeof(metadata(stack)) <: Dict
        @test refdims(stack) == ()

        # Loads child as a regular GeoArray
        @test typeof(stack[:albedo]) <: GeoArray{Union{Missing,Float32},3} 
        @test typeof(stack[:albedo, 2, 3, 1]) <: Float32
        @test typeof(stack[:albedo, :, 3, 1]) <: GeoArray{Union{Missing,Float32},1}
        @test typeof(dims(stack, :albedo)) <: Tuple{<:Lon,<:Lat,<:Time}
        @test typeof(keys(stack)) == NTuple{131,Symbol}
        @test first(keys(stack)) == :abso4
        @test typeof(metadata(stack, :albedo)) <: Dict
        @test metadata(stack, :albedo)["institution"] == "Max-Planck-Institute for Meteorology"

        # Test some DimensionalData.jl tools work
        # Time dim should be reduced to length 1 by mean
        @test axes(mean(stack[:albedo, Lat(1:20)], dims=Time)) == (Base.OneTo(192), Base.OneTo(20), Base.OneTo(1))
        array = stack[:albedo][Time(4:6), Lon(1), Lat(2)] 
        @test array == stack[:albedo, Time(4:6), Lon(1), Lat(2)] 
        size(array) == (3,)
    end

    @testset "copy" begin
        array = GeoArray(stack[:albedo])
        copy!(array, stack, :albedo)
    end

    @testset "indexing" begin
        ncmultistack = NCstack([geturl(ncsingle)])
        ncmultistack = NCstack((geturl(ncsingle),))
        @test typeof(ncmultistack[:tos]) <: GeoArray{Union{Missing,Float32},3}
        @test typeof(ncmultistack[:tos, Time(1)]) <: GeoArray{Union{Missing,Float32},2}
        @test typeof(ncmultistack[:tos, Lat(1), Time(1)]) <: GeoArray{Union{Missing,Float32},1}
        @test typeof(ncmultistack[:tos, 8, 30, 10]) <: Float32
    end

    @testset "conversion to GeoStack" begin
        stack = GeoStack(stack)
        @test Symbol.(Tuple(keys(stack))) == keys(stack)
        smallstack = GeoStack(stack; keys=(:albedo, :evap, :runoff))
        @test keys(smallstack) == (:albedo, :evap, :runoff)
    end
end


