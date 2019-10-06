path = geturl("https://download.osgeo.org/geotiff/samples/gdal_eg/cea.tif")
path = geturl("https://download.osgeo.org/geotiff/samples/usgs/c41078a1.tif")
path = geturl("https://download.osgeo.org/geotiff/samples/usgs/i30dem.tif")

@testset "array" begin
    array = GDALarray(path)

    @testset "array properties" begin
        @test size(array) == (514, 515, 1)
        @test typeof(array) <: GDALarray{UInt8,3}
    end

    @testset "dimensions" begin
        @test length(val(dims(dims(array), Lon))) == 514
        @test ndims(array) == 3
        @test typeof(dims(array)) <: Tuple{<:Lon,<:Lat,<:Band}
        @test refdims(array) == ()
        # @test bounds(array) 
    end

    @testset "other fields" begin
        @test window(array) == ()
        @test missingval(array) == -1.0e10
        @test typeof(metadata(array)) <: NamedTuple
        @test metadata(array).filepath == "cea.tif"
        @test name(array) == Symbol("")
    end

    @testset "indexing" begin 
        @test typeof(array[Band(1)]) <: GeoArray{UInt8,2} 
        @test typeof(array[Lat(1), Band(1)]) <: GeoArray{UInt8,1} 
        @test typeof(array[Lon(1), Band(1)]) <: GeoArray{UInt8,1}
        # Doesn't handle returning a single value
        @test_broken typeof(array[Lon(1), Lat(1), Band(1)]) <: UInt8 
        @test_broken typeof(array[1, 1, 1]) <: UInt8
    end

    @testset "selectors" begin
        a = array[Lon(At(3)), Lat(:), Band(1)]
        @test typeof(a) <: GeoArray{UInt8,2}
        @test bounds(a) == ((20, 30), (20, 30))
        # Doesn't handle returning a single value
        # a = array[Lon(At(20), Lat(Near(10), Band(1)]) <: UInt8
    end

    @testset "conversion to GeoArray" begin
        geoarray = array[Lon(1:50), Lat(1:1), Band(1)]
        @test size(geoarray) == (50, 1)
        @test eltype(geoarray) <: UInt8
        @time typeof(geoarray) <: GeoArray{UInt8,1} 
        @test typeof(dims(geoarray)) <: Tuple{<:Lon,Lat}
        @test typeof(refdims(geoarray)) <: Tuple{<:Band} 
        @test metadata(geoarray) == metadata(array)
        @test missingval(geoarray) == -1.0e10
        @test name(array) == Symbol("")
    end
end

@testset "stack" begin
    gdalstack = GDALstack((a=geturl(gdal_url),))

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
        array = zero(GeoArray(gdalstack[:a]))
        copy!(array, gdalstack, :a)
    end
end

"Driver: GTiff/GeoTIFF
Files: /home/raf/CESAR/Raster/limited_growth/limited_growth_2016_01.tif
Size is 3856, 1624
Coordinate System is:
PROJCS[\"unnamed\",
    GEOGCS[\"WGS 84\",
        DATUM[\"WGS_1984\",
            SPHEROID[\"WGS 84\",6378137,298.257223563,
                AUTHORITY[\"EPSG\",\"7030\"]],
            AUTHORITY[\"EPSG\",\"6326\"]],
        PRIMEM[\"Greenwich\",0],
        UNIT[\"degree\",0.0174532925199433],
        AUTHORITY[\"EPSG\",\"4326\"]],
    PROJECTION[\"Cylindrical_Equal_Area\"],
    PARAMETER[\"standard_parallel_1\",30],
    PARAMETER[\"central_meridian\",0],
    PARAMETER[\"false_easting\",0],
    PARAMETER[\"false_northing\",0],
    UNIT[\"metre\",1,
        AUT
HORITY[\"EPSG\",\"9001\"]]]
Origin = (-17367530.445161368697882,7314540.795860165730119)
Pixel Size = (9008.055210145937963,-9008.055167315476865)
Metadata:
  AREA_OR_POINT=Area
Image Structure Metadata:
  COMPRESSION=LZW
  INTERLEAVE=BAND
Corner Coordinates:
Upper Left  (-17367530.445, 7314540.796) (180d 0' 0.00\"W, 85d 2'40.43\"N)
Lower Left  (-1
7367530.445,-7314540.796) (180d 0' 0.00\"W, 85d 2'40.43\"S)
Upper Right (17367530.445, 7314540.796) (180d 0' 0.00\"E, 8
5d 2'40.43\"N)
Lower Right (17367530.445,-7314540.796) (180d 0' 0.00\"E, 85d 2'40.43\"S)
Center      (   0.0000000,  -
0.0000000) (  0d 0' 0.01\"E,  0d 0' 0.00\"S)
Band 1 Block=3856x1 Type=Float32, ColorInterp=Gray
  Min=-59.869 Max=11.4
85 
  Minimum=-59.869, Maximum=11.485, Mean=-13.064, StdDev=13.656
  NoData Value=-3.39999999999999996e+38
  Metadata
:
    STATISTICS_MAXIMUM=11.484945361028
    STATISTICS_MEAN=-13.064166792583
    STATISTICS_MINIMUM=-59.868704210423

    STATISTICS_STDDEV=13.656090472918
"

path = geturl("https://download.osgeo.org/geotiff/samples/usgs/o41078a7.tif")

path = "/home/raf/CESAR/Raster/limited_growth/limited_growth_2016_01.tif"

using GeoArrays
using Plots
geoarray = GeoArrays.read(path)
geoarray.f
geoarray.crs
coords(geoarray, [1,1]) ./ 9008.055210145937963
coords(geoarray, [size(geoarray)[1:2]...]) .* 1e-5
indices(geoarray, [440720.0, 3.75132e6])
coords(geoarray)
using ArchGDAL
ArchGDAL.registerdrivers() do 
    ArchGDAL.read(path) do ds
        ArchGDAL.gdalinfo(ds)
    end
end

    bounds(array)
    map((s, b)-> (b[2]-b[1]) / s, size(array), bounds(array))
    size(array)
    3856 * 0.09008055210145939
    360/3856
    missingval(array)
    display(string(crs(array)))
    plotly()
    plot(array)
