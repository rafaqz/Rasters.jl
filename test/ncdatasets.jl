ncexamples = "https://www.unidata.ucar.edu/software/netcdf/examples/"
ncsingle = geturl(joinpath(ncexamples, "tos_O1_2001-2002.nc"))
ncmulti = geturl(joinpath(ncexamples, "test_echam_spectral.nc"))

@testset "NCarray" begin
    ncarray = NCDarray(ncsingle)

    @testset "array properties" begin
        @test size(ncarray) == (180, 170, 24)
        @test typeof(ncarray) <: NCDarray{Union{Missing,Float32},3}
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
        @test typeof(metadata(ncarray)) <: Dict # TODO make this a namedtuple
        @test_broken metadata(ncarray).filepath == "tos_O1_2001-2002.nc"
        @test name(ncarray) == :tos
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
        @test name(geoarray) == :tos
    end

end

@testset "NCDstack" begin
    ncstack = NCDstack(geturl(ncmulti))

    @testset "load ncstack" begin
        @test typeof(ncstack) <: NCDstack{String}
        @test ismissing(missingval(ncstack))
        @test typeof(metadata(ncstack)) <: Dict
        @test refdims(ncstack) == ()
        # Loads child as a regular GeoArray
        @test typeof(ncstack[:albedo]) <: GeoArray{Union{Missing,Float32},3} 
        @test typeof(ncstack[:albedo, 2, 3, 1]) <: Float32
        @test typeof(ncstack[:albedo, :, 3, 1]) <: GeoArray{Union{Missing,Float32},1}
        @test typeof(dims(ncstack, :albedo)) <: Tuple{<:Lon,<:Lat,<:Time}
        @test typeof(keys(ncstack)) == NTuple{131,Symbol}
        @test first(keys(ncstack)) == :abso4
        @test typeof(metadata(ncstack, :albedo)) <: Dict
        @test metadata(ncstack, :albedo)["institution"] == "Max-Planck-Institute for Meteorology"
        # Test some DimensionalData.jl tools work
        # Time dim should be reduced to length 1 by mean
        @test axes(mean(ncstack[:albedo, Lat(1:20)], dims=Time)) == (Base.OneTo(192), Base.OneTo(20), Base.OneTo(1))
        geoarray = ncstack[:albedo][Time(4:6), Lon(1), Lat(2)] 
        @test geoarray == ncstack[:albedo, Time(4:6), Lon(1), Lat(2)] 
        size(geoarray) == (3,)
    end

    @testset "copy" begin
        geoarray = GeoArray(ncstack[:albedo])
        copy!(geoarray, ncstack, :albedo)
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

    @testset "conversion to GeoStack" begin
        geostack = GeoStack(ncstack)
        @test Symbol.(Tuple(keys(geostack))) == keys(ncstack)
        smallstack = GeoStack(ncstack; keys=(:albedo, :evap, :runoff))
        @test keys(smallstack) == (:albedo, :evap, :runoff)
    end

end


