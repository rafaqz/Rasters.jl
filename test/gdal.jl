using ArchGDAL, GeoData, Test, Statistics, Dates
using GeoData: window
include("test_utils.jl")

path = geturl("https://download.osgeo.org/geotiff/samples/gdal_eg/cea.tif")

@testset "array" begin
    gdalarray = GDALarray(path)

    @testset "array properties" begin
        @test size(gdalarray) == (514, 515, 1)
        @test gdalarray isa GDALarray{UInt8,3}
    end

    @testset "dimensions" begin
        @test length(val(dims(dims(gdalarray), Lon))) == 514
        @test ndims(gdalarray) == 3
        @test dims(gdalarray) isa Tuple{<:Lon,<:Lat,<:Band}
        @test refdims(gdalarray) == ()
        @test_broken bounds(gdalarray) 
    end

    @testset "other fields" begin
        @test window(gdalarray) == ()
        @test missingval(gdalarray) == -1.0e10
        @test metadata(gdalarray) isa GDALmetadata
        @test basename(metadata(gdalarray).val["filepath"]) == "cea.tif"
        @test name(gdalarray) == "Unnamed"
    end

    @testset "indexing" begin 
        @test gdalarray[Band(1)] isa GeoArray{UInt8,2} 
        @test gdalarray[Lat(1), Band(1)] isa GeoArray{UInt8,1} 
        @test gdalarray[Lon(1), Band(1)] isa GeoArray{UInt8,1}
        @test gdalarray[Lon(1), Lat(1), Band(1)] isa UInt8 
        @test gdalarray[1, 1, 1] isa UInt8
    end

    @testset "selectors" begin
        geoarray = gdalarray[Lat(Near(3)), Lon(:), Band(1)]
        @test geoarray isa GeoArray{UInt8,1}
        @test gdalarray[Lon(10), Lat(10), Band(1)] == 0x73
    end

    @testset "conversion to GeoArray" begin
        geoarray = gdalarray[Lon(1:50), Lat(1:1), Band(1)]
        @test size(geoarray) == (50, 1)
        @test eltype(geoarray) <: UInt8
        @time geoarray isa GeoArray{UInt8,1} 
        @test dims(geoarray) isa Tuple{<:Lon,Lat}
        @test refdims(geoarray) isa Tuple{<:Band} 
        @test metadata(geoarray) == metadata(gdalarray)
        @test missingval(geoarray) == -1.0e10
        @test name(geoarray) == "Unnamed"
    end

    # Works but saved raster has no geotransform so can't be loaded
    # @testset "save" begin
    #     geoarray = gdalarray[Band(1)]
    #     filename = tempname()
    #     GeoData.write(filename, GDALarray, geoarray)
    #     saved = GeoArray(GDALarray(filename))
    #     # 1 bands is added again on save
    #     @test size(saved) != size(geoarray)
    #     @test size(saved[Band(1)]) == size(geoarray)
    #     @test refdims(saved) == refdims(geoarray)
    #     @test missingval(saved) === missingval(geoarray)
    #     @test metadata(saved) == metadata(geoarray)
    #     @test GeoData.name(saved) == GeoData.name(geoarray)
    #     @test all(metadata.(dims(saved)) .== metadata.(dims(geoarray)))
    #     @test all(DimensionalData.grid.(dims(saved)) .== DimensionalData.grid.(dims(geoarray)))
    #     @test typeof(dims(saved)) == typeof(dims(geoarray))
    #     @test_broken val(dims(saved)[3]) == val(dims(geoarray)[3])
    #     @test_broken all(val.(dims(saved)) .== val.(dims(geoarray)))
    #     @test all(metadata.(dims(saved)) .== metadata.(dims(geoarray)))
    #     @test all(data(saved) .=== data(geoarray))
    #     @test typeof(saved) == typeof(geoarray)
    #     geoarray = gdalarray
    #     GeoData.write(filename, GDALarray, geoarray)
    #     saved = GeoArray(GDALarray(filename))
    #     @test size(saved) == size(geoarray)
    # end

end

@testset "stack" begin
    gdalstack = GDALstack((a=path, b=path));

    @testset "child array properties" begin
        @test size(gdalstack[:a]) == (514, 515, 1)
        @test gdalstack[:a] isa GDALarray{UInt8,3}
    end

    @testset "indexing" begin
        @test gdalstack[:a][Lat(2:3), Lon(1), Band(1)] == [0x00, 0x6b]
        @test gdalstack[:a][Lat(1), Lon(1), Band(1)] == 0x00
    end

    @testset "window" begin
        windowedstack = GDALstack((a=path, b=path); window=(Lat(1:5), Lon(1:5), Band(1)))
        @test window(windowedstack) == (Lat(1:5), Lon(1:5), Band(1))
        windowedarray = GeoArray(windowedstack[:a])
        @test windowedarray isa GeoArray{UInt8,2}
        @test length.(dims(windowedarray)) == (5, 5)
        @test size(windowedarray) == (5, 5)
        @test windowedarray[1:3, 2:2] == reshape([0x00, 0x00, 0x00], 3, 1)
        @test windowedarray[1:3, 2] == [0x00, 0x00, 0x00]
        @test windowedarray[1, 2] == 0x00
        windowedstack = GDALstack((a=path, b=path); window=(Lat(1:5), Lon(1:5), Band(1:1)))
        windowedarray = windowedstack[:b]
        @test windowedarray[1:3, 2:2, 1:1] == reshape([0x00, 0x00, 0x00], 3, 1, 1)
        @test windowedarray[1:3, 2:2, 1] == reshape([0x00, 0x00, 0x00], 3, 1)
        @test windowedarray[1:3, 2, 1] == [0x00, 0x00, 0x00]
        @test windowedarray[1, 2, 1] == 0x00
        windowedstack = GDALstack((a=path, b=path); window=Band(1))
        windowedarray = GeoArray(windowedstack[:b])
        @test windowedarray[1:3, 2:2] == reshape([0x00, 0x00, 0x00], 3, 1)
        @test windowedarray[1:3, 2] == [0x00, 0x00, 0x00]
        @test windowedarray[1, 2] == 0x00
    end


    # Stack Constructors
    @testset "conversion to GeoStack" begin
        stack = GeoStack(gdalstack)
        @test Symbol.(Tuple(keys(gdalstack))) == keys(stack)
        smallstack = GeoStack(gdalstack; keys=(:a,))
        keys(smallstack) == (:a,)
    end

    if VERSION > v"1.1-"
        @testset "copy" begin
            geoarray = zero(GeoArray(gdalstack[:a]))
            copy!(geoarray, gdalstack, :a)
            # First wrap with GeoArray() here or == loads from disk for each cell.
            # we need a general way of avoiding this in all disk-based sources
            @test geoarray == GeoArray(gdalstack[:a])
        end
    end

    # @testset "save" begin
    #     geoarray = GeoArray(geostack[:a])
    #     filename = tempname()
    #     write(filename, GDALarray, gdalstack)
    #     base, ext = splitext(filename)
    #     filename_b = string(base, "_b", ext)
    #     saved = GeoArray(GrdArray(filename_b))
    #     @test saved == geoarray
    # end

end
