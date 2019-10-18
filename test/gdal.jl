path = geturl("https://download.osgeo.org/geotiff/samples/usgs/c41078a1.tif")
path = geturl("https://download.osgeo.org/geotiff/samples/usgs/i30dem.tif")
path = geturl("https://download.osgeo.orgtgeotiff/samples/gdal_eg/cea.tif")

@testset "array" begin
    A = GDALarray(path)

    @testset "array properties" begin
        @test size(A) == (514, 515, 1)
        @test typeof(A) <: GDALarray{UInt8,3}
    end

    @testset "dimensions" begin
        @test length(val(dims(dims(A), Lon))) == 514
        @test ndims(A) == 3
        @test typeof(dims(A)) <: Tuple{<:Lon,<:Lat,<:Band}
        @test refdims(A) == ()
        @test_broken bounds(A) 
    end

    @testset "other fields" begin
        @test window(A) == ()
        @test missingval(A) == -1.0e10
        @test typeof(metadata(A)) <: NamedTuple
        @test metadata(A).filepath == "cea.tif"
        @test name(A) == Symbol("")
    end

    @testset "indexing" begin 
        @test typeof(A[Band(1)]) <: GeoArray{UInt8,2} 
        @test typeof(A[Lat(1), Band(1)]) <: GeoArray{UInt8,1} 
        @test typeof(A[Lon(1), Band(1)]) <: GeoArray{UInt8,1}
        # Doesn't handle returning a single value
        @test_broken typeof(A[Lon(1), Lat(1), Band(1)]) <: UInt8 
        @test_broken typeof(A[1, 1, 1]) <: UInt8
    end

    @testset "selectors" begin
        a = A[Lat(Near(3)), Lon(:), Band(1)]
        @test typeof(a) <: GeoArray{UInt8,1}
        # @test bounds(a) == ()
        # Doesn't handle returning a single value
        # a = A[Lon(At(20), Lat(Near(10), Band(1)]) <: UInt8
    end

    @testset "conversion to GeoArray" begin
        geoarray = A[Lon(1:50), Lat(1:1), Band(1)]
        @test size(geoarray) == (50, 1)
        @test eltype(geoarray) <: UInt8
        @time typeof(geoarray) <: GeoArray{UInt8,1} 
        @test typeof(dims(geoarray)) <: Tuple{<:Lon,Lat}
        @test typeof(refdims(geoarray)) <: Tuple{<:Band} 
        @test metadata(geoarray) == metadata(A)
        @test missingval(geoarray) == -1.0e10
        @test name(geoarray) == Symbol("")
    end

end

@testset "stack" begin
    gdalstack = GDALstack((a=path,))

    # Broken: ArchGDAL read() indexing is non-standard
    # band is first when it is last in the returned array.
    # At least one unitrange is required for dispatch.
    @test_broken gdalstack[:a][Lat([2,3]), Lon(1), Band(1)]
    @test_broken gdalstack[:a][Lat(1), Lon(1), Band(1)]

    # Stack Constructors
    @testset "conversion to GeoStack" begin
        stack = GeoStack(gdalstack)
        @test Symbol.(Tuple(keys(gdalstack))) == keys(stack)
        smallstack = GeoStack(gdalstack; keys=(:a,))
        keys(smallstack) == (:a,)
    end

    @testset "copy" begin
        A = zero(GeoArray(gdalstack[:a]))
        copy!(A, gdalstack, :a)
        # First wrap with GeoArray() here or == loads from disk for each cell.
        # we need a general way of avoiding this in all disk-based sources
        @test A == GeoArray(gdalstack[:a])
    end
end
