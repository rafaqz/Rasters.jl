using ArchGDAL, GeoData, Test, Statistics, Dates, Plots, NCDatasets
using GeoData: window, mode, span, sampling, name, dims2indices

include(joinpath(dirname(pathof(GeoData)), "../test/test_utils.jl"))

path = geturl("https://download.osgeo.org/geotiff/samples/gdal_eg/cea.tif")

@testset "GDALarray" begin
    gdalarray = GDALarray(path; usercrs=EPSG(4326), name="test");

    @testset "array properties" begin
        @test size(gdalarray) == (514, 515, 1)
        @test gdalarray isa GDALarray{UInt8,3}
    end

    coord = first.(bounds(dims(gdalarray, (Lon, Lat))))
    ArchGDAL.reproject([20, 20], EPSG(4326), ProjString("+proj=cea"); order=:trad)

    @testset "dimensions" begin
        @test length(val(dims(dims(gdalarray), Lon))) == 514
        @test ndims(gdalarray) == 3
        @test dims(gdalarray) isa Tuple{<:Lon,<:Lat,<:Band}
        @test mode(dims(gdalarray, Band)) == DimensionalData.Categorical(Ordered())
        @test span.(mode.(dims(gdalarray, (Lat, Lon)))) == 
            (Regular(-60.02213698319351), Regular(60.02213698319374))
        @test sampling.(mode.(dims(gdalarray, (Lat, Lon)))) == 
            (Intervals(Start()), Intervals(Start()))
        @test refdims(gdalarray) == ()
        @test bounds(gdalarray) == ((-28493.166784412522, 2358.2116249490587), 
                                    (4.22503316539283e6, 4.255944565939175e6), 
                                    (1, 1))
    end

    @testset "other fields" begin
        @test missingval(gdalarray) == -1.0e10
        @test metadata(gdalarray) isa GDALarrayMetadata
        @test basename(metadata(gdalarray).val["filepath"]) == "cea.tif"
        @test name(gdalarray) == "test"
        @test label(gdalarray) == "test"
        @test units(gdalarray) == nothing
        @test usercrs(dims(gdalarray, Lat)) == EPSG(4326)
        @test usercrs(dims(gdalarray, Lon)) == EPSG(4326)
        @test usercrs(gdalarray) == EPSG(4326)
        @test usercrs(gdalarray[Lat(1)]) == EPSG(4326)
        @test_throws ErrorException usercrs(gdalarray[Lat(1), Lon(1)])
        wkt = WellKnownText(GeoFormatTypes.CRS(), 
          "PROJCS[\"unnamed\",GEOGCS[\"NAD27\",DATUM[\"North_American_Datum_1927\",SPHEROID[\"Clarke 1866\",6378206.4,294.978698213898,AUTHORITY[\"EPSG\",\"7008\"]],AUTHORITY[\"EPSG\",\"6267\"]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4267\"]],PROJECTION[\"Cylindrical_Equal_Area\"],PARAMETER[\"standard_parallel_1\",33.75],PARAMETER[\"central_meridian\",-117.333333333333],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH]]")
        @test crs(dims(gdalarray, Lat)) == wkt
        @test crs(dims(gdalarray, Lon)) == wkt
        @test crs(gdalarray) == wkt
        @test crs(gdalarray[Lat(1)]) == wkt
        @test_throws ErrorException crs(gdalarray[Lat(1), Lon(1)])
    end

    @testset "indexing" begin 
        @test gdalarray[Band(1)] isa GeoArray{UInt8,2} 
        @test gdalarray[Lat(1), Band(1)] isa GeoArray{UInt8,1} 
        @test gdalarray[Lon(1), Band(1)] isa GeoArray{UInt8,1}
        @test gdalarray[Lon(1), Lat(1), Band(1)] isa UInt8 
        @test gdalarray[1, 1, 1] isa UInt8
        @test gdalarray[Lat(Contains(33.92)), Lon(Contains(-117.56)), Band(1)] == parent(gdalarray)[125, 45, 1] == 0xbd
    end

    @testset "methods" begin 
        mean(gdalarray; dims=Lat) == mean(data(gdalarray); dims=2)
    end

    @testset "selectors" begin
        # TODO verify the value with R/gdal etc
        @test gdalarray[Lat(Contains(33.8)), Lon(Contains(-117.5)), Band(1)] isa UInt8
        @test gdalarray[Lat(Between(33.7, 33.9)), Band(1)] isa GeoArray
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
        @test name(geoarray) == "test"
    end

    @testset "save" begin
        gdalarray = GDALarray(path; usercrs=EPSG(4326));

        @testset "2d" begin
            geoarray1 = gdalarray[Band(1)]
            filename = tempname() * ".tif"
            write(filename, GDALarray, geoarray1)
            saved1 = GDALarray(filename; usercrs=EPSG(4326))[Band(1)];
            dims(saved1)
            @test saved1 == geoarray1
            @test typeof(saved1) == typeof(geoarray1)
            @test val(dims(saved1, Lon)) == val(dims(geoarray1, Lon))
            @test val(dims(saved1, Lat)) == val(dims(geoarray1, Lat))
            @test all(metadata.(dims(saved1)) .== metadata.(dims(geoarray1)))
            @test missingval(saved1) === missingval(geoarray1)
            @test refdims(saved1) == refdims(geoarray1)
        end
        
        @testset "3d, with subsetting" begin
            geoarray2 = gdalarray[Lat(Between(33.7, 33.9)), 
                                  Lon(Between(-117.6, -117.4))]
            filename2 = tempname() * ".img"
            write(filename2, GDALarray, geoarray2)
            saved2 = GeoArray(GDALarray(filename2; usercrs=EPSG(4326)))
            @test size(saved2) == size(geoarray2) == length.(dims(saved2)) == length.(dims(geoarray2))
            @test refdims(saved2) == refdims(geoarray2)
            #TODO test a file with more metadata
            @test metadata(saved2)["filepath"] == filename2
            @test missingval(saved2) === missingval(geoarray2)
            @test GeoData.name(saved2) == GeoData.name(geoarray2)
            @test step(mode(dims(saved2, Lat))) ≈ step(mode(dims(geoarray2, Lat)))
            @test step(mode(dims(saved2, Lon))) ≈ step(mode(dims(geoarray2, Lon)))
            @test typeof(dims(saved2)) == typeof(dims(geoarray2))
            @test all(val(dims(saved2, Band)) .≈ val(dims(geoarray2, Band)))
            @test all(val(dims(saved2, Lon)) .≈ val(dims(geoarray2, Lon)))
            @test all(val(dims(saved2, Lat)) .≈ val(dims(geoarray2, Lat)))
            @test all(metadata.(dims(saved2)) .== metadata.(dims(geoarray2)))
            @test data(saved2) == data(geoarray2)
            @test typeof(saved2) == typeof(geoarray2)
            filename3 = tempname() * ".img"
            geoarray3 = cat(gdalarray[Band(1)], gdalarray[Band(1)], gdalarray[Band(1)]; dims=Band(1:3))
            write(filename3, GDALarray, geoarray3)
            saved3 = GeoArray(GDALarray(filename3; usercrs=EPSG(4326)))
            @test saved3 == geoarray3
            @test val(dims(saved3, Band)) == 1:3
        end

        @testset "resave current" begin
            filename = tempname() * ".tiff"
            write(filename, gdalarray)
            gdalarray2 = GDALarray(filename)
            write(gdalarray2)
            @test convert(GeoArray, GDALarray(filename)) == convert(GeoArray, gdalarray2)
        end

        @testset "to grd" begin
            write("testgrd", GRDarray, gdalarray)
            grdarray = GRDarray("testgrd")
            @test crs(grdarray) == convert(ProjString, crs(gdalarray))
            @test bounds(grdarray) == (bounds(gdalarray))
            @test val(dims(grdarray, Lat)) == reverse(val(dims(gdalarray, Lat)))
            @test val(dims(grdarray, Lon)) ≈ val(dims(gdalarray, Lon))
            @test GeoArray(grdarray) == GeoArray(gdalarray)
            @test bounds(grdarray) == bounds(gdalarray)
        end

        @testset "to netcdf" begin
            filename2 = tempname()
            write(filename2, NCDarray, gdalarray[Band(1)])
            saved = GeoArray(NCDarray(filename2))
            @test size(saved) == size(gdalarray[Band(1)])
            @test saved ≈ reverse(gdalarray[Band(1)]; dims=Lat)
            @test val(dims(saved, Lon)) ≈ val(dims(gdalarray, Lon)) .+ 0.5step(dims(saved, Lon))
            @test val(dims(saved, Lat)) ≈ reverse(val(dims(gdalarray, Lat))) .+ 0.5step(dims(saved, Lat))
            @test all(map(isapprox, bounds(saved, Lon), bounds(gdalarray, Lon)))
            @test all(map(isapprox, bounds(saved, Lat), bounds(gdalarray, Lat)))
        end

    end

    @testset "plot" begin
        # TODO write some tests for this
        gdalarray |> plot
    end

end

@testset "GDAL stack" begin
    gdalstack = GDALstack((a=path, b=path));

    @testset "child array properties" begin
        @test size(gdalstack[:a]) == (514, 515, 1)
        @test gdalstack[:a] isa GeoArray{UInt8,3}
    end

    @testset "indexing" begin
        @test gdalstack[:a][Lat(2:3), Lon(1), Band(1)] == [0x00, 0x6b]
        @test gdalstack[:a][Lat(1), Lon(1), Band(1)] == 0x00
        @test gdalstack[:b, Band(1)] == gdalstack[:b][Band(1)]
        @test typeof(gdalstack[:b, Band(1)]) == typeof(gdalstack[:b][Band(1)])
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
        windowedstack = GDALstack((a=path, b=path); window=(Band(1),))
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
        @test keys(smallstack) == (:a,)
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

    @testset "save" begin
        geoarray = GeoArray(gdalstack[:a])
        filename = tempname()
        write(filename, GDALarray, gdalstack)
        base, ext = splitext(filename)
        filename_b = string(base, "_b", ext)
        saved = GeoArray(GDALarray(filename_b))
        @test saved == geoarray
    end

end

@testset "GDAL series" begin
    series = GeoSeries([path, path], (Ti,); childtype=GDALarray, childkwargs=(usercrs=EPSG(4326), name="test"))
    @test GeoArray(series[Ti(1)]) == GeoArray(GDALarray(path; usercrs=EPSG(4326), name="test"))

    gdalstack = GDALstack((a=path, b=path); childtype=GDALarray, childkwargs=(usercrs=EPSG(4326),))
    series = GeoSeries([gdalstack, gdalstack], (Ti,))
    @test series[1].childkwargs == gdalstack.childkwargs
    # Rebuild the series by wrapping the GDALarray data in Array.
    # `modify` forces `rebuild` on all containers as in-Memory variants
    modified_series = modify(Array, series)
    @test typeof(modified_series) <: GeoSeries{<:GeoStack{<:NamedTuple{(:a,:b),<:Tuple{<:GeoArray{UInt8,3,<:Tuple,<:Tuple,<:Array{UInt8,3}},Vararg}}}}
end

nothing
