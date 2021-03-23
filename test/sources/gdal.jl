using GeoData, Test, Statistics, Dates, Plots, DimensionalData
import ArchGDAL, NCDatasets
using GeoData: window, mode, span, sampling, name, bounds

include(joinpath(dirname(pathof(GeoData)), "../test/test_utils.jl"))

path = maybedownload("https://download.osgeo.org/geotiff/samples/gdal_eg/cea.tif")

@testset "GDALarray" begin
    gdalarray = geoarray(path; mappedcrs=EPSG(4326), name=:test);

    @testset "open" begin
        @test open(A -> A[Y=1], gdalarray) == gdalarray[:, 1, :]
    end

    @testset "array properties" begin
        @test size(gdalarray) == (514, 515, 1)
        @test gdalarray isa GDALarray{UInt8,3}
    end

    @testset "dimensions" begin
        @test length(dims(gdalarray, X)) == 514
        @test ndims(gdalarray) == 3
        @test dims(gdalarray) isa Tuple{<:X,<:Y,<:Band}
        @test mode(gdalarray, Band) == DimensionalData.Categorical(Ordered())
        @test span(gdalarray, (Y, X)) == 
            (Regular(-60.02213698319351), Regular(60.02213698319374))
        @test sampling(gdalarray, (Y, X)) == 
            (Intervals(Start()), Intervals(Start()))
        @test refdims(gdalarray) == ()
        # Bounds calculated in python using rasterio
        @test all(bounds(gdalarray, Y) .≈ (4224973.143255847, 4255884.5438021915))
        @test all(bounds(gdalarray, X) .≈ (-28493.166784412522, 2358.211624949061))
    end

    @testset "other fields" begin
        # This file has an inorrect missing value
        @test missingval(gdalarray) == nothing
        @test metadata(gdalarray) isa Metadata{:GDAL}
        @test basename(metadata(gdalarray).val[:filepath]) == "cea.tif"
        @test name(gdalarray) == :test
        @test label(gdalarray) == "test"
        @test units(gdalarray) == nothing
        @test mappedcrs(dims(gdalarray, Y)) == EPSG(4326)
        @test mappedcrs(dims(gdalarray, X)) == EPSG(4326)
        @test mappedcrs(gdalarray) == EPSG(4326)
        @test mappedcrs(gdalarray[Y(1)]) == EPSG(4326)
        @test_throws ErrorException mappedcrs(gdalarray[Y(1), X(1)])
        wkt = WellKnownText(GeoFormatTypes.CRS(), 
          "PROJCS[\"unnamed\",GEOGCS[\"NAD27\",DATUM[\"North_American_Datum_1927\",SPHEROID[\"Clarke 1866\",6378206.4,294.978698213898,AUTHORITY[\"EPSG\",\"7008\"]],AUTHORITY[\"EPSG\",\"6267\"]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4267\"]],PROJECTION[\"Cylindrical_Equal_Area\"],PARAMETER[\"standard_parallel_1\",33.75],PARAMETER[\"central_meridian\",-117.333333333333],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH]]")
        @test crs(dims(gdalarray, Y)) == wkt
        @test crs(dims(gdalarray, X)) == wkt
        @test crs(gdalarray) == wkt
        @test crs(gdalarray[Y(1)]) == wkt
        @test_throws ErrorException crs(gdalarray[Y(1), X(1)])
    end

    @testset "indexing" begin 
        @test gdalarray[Band(1)] isa GeoArray{UInt8,2} 
        @test gdalarray[Y(1), Band(1)] isa GeoArray{UInt8,1} 
        @test gdalarray[X(1), Band(1)] isa GeoArray{UInt8,1}
        @test gdalarray[X(1), Y(1), Band(1)] isa UInt8 
        @test gdalarray[1, 1, 1] isa UInt8
        # Value indexed in python with rasterio
        @test gdalarray[Y(21), X(21), Band(1)] == 82 
    end

    @testset "methods" begin 
        mean(gdalarray; dims=Y) == mean(data(gdalarray); dims=2)
    end

    @testset "selectors" begin
        # TODO verify the value with R/gdal etc
        @test gdalarray[Y(Contains(33.8)), X(Contains(-117.5)), Band(1)] isa UInt8
        @test gdalarray[Y(Between(33.7, 33.9)), Band(1)] isa GeoArray
    end

    @testset "conversion to GeoArray" begin
        geoA = gdalarray[X(1:50), Y(1:1), Band(1)]
        @test size(geoA) == (50, 1)
        @test eltype(geoA) <: UInt8
        @time geoA isa GeoArray{UInt8,1} 
        @test dims(geoA) isa Tuple{<:X,Y}
        @test refdims(geoA) isa Tuple{<:Band} 
        @test metadata(geoA) == metadata(gdalarray)
        @test missingval(geoA) == nothing
        @test name(geoA) == :test
    end

    @testset "save" begin
        gdalarray = GDALarray(path; mappedcrs=EPSG(4326), name=:test);

        @testset "2d" begin
            geoA = gdalarray[Band(1)]
            filename = tempname() * ".asc"
            write(filename, geoA)
            saved1 = GDALarray(filename; mappedcrs=EPSG(4326))[Band(1)];
            @test all(saved1 .== saved1 ≈ geoA)
            @test typeof(saved1) !== typeof(geoA)
            @test val(dims(saved1, X)) ≈ val(dims(geoA, X))
            @test val(dims(saved1, Y)) ≈ val(dims(geoA, Y))
            @test all(metadata.(dims(saved1)) .== metadata.(dims(geoA)))
            @test metadata(dims(saved1)[1]) == metadata(dims(geoA)[1])
            @test missingval(saved1) === missingval(geoA) 
            @test refdims(saved1) == refdims(geoA) 
        end
        
        @testset "3d, with subsetting" begin
            geoA2 = gdalarray[Y(Between(33.7, 33.9)), 
                                  X(Between(-117.6, -117.4))]
            filename2 = tempname() * ".tif"
            write(filename2, geoA2)
            saved2 = GeoArray(GDALarray(filename2; name=:test, mappedcrs=EPSG(4326)))
            @test size(saved2) == size(geoA2) == length.(dims(saved2)) == length.(dims(geoA2))
            @test refdims(saved2) == refdims(geoA2)
            #TODO test a file with more metadata
            @test val(metadata(saved2))[:filepath] == filename2
            @test missingval(saved2) === missingval(geoA2)
            @test GeoData.name(saved2) == GeoData.name(geoA2)
            @test step(mode(dims(saved2, Y))) ≈ step(mode(dims(geoA2, Y)))
            @test step(mode(dims(saved2, X))) ≈ step(mode(dims(geoA2, X)))
            @test typeof(dims(saved2)) == typeof(dims(geoA2))
            @test all(val(dims(saved2, Band)) .≈ val(dims(geoA2, Band)))
            @test all(val(dims(saved2, X)) .≈ val(dims(geoA2, X)))
            @test all(val(dims(saved2, Y)) .≈ val(dims(geoA2, Y)))
            @test all(metadata.(dims(saved2)) .== metadata.(dims(geoA2)))
            @test data(saved2) == data(geoA2)
            @test typeof(saved2) == typeof(geoA2)
            filename3 = tempname() * ".tif"
            geoA3 = cat(gdalarray[Band(1)], gdalarray[Band(1)], gdalarray[Band(1)]; dims=Band(1:3))
            write(filename3, geoA3)
            saved3 = GeoArray(GDALarray(filename3; mappedcrs=EPSG(4326)))
            @test all(saved3 .== geoA3)
            @test val(dims(saved3, Band)) == 1:3
        end

        @testset "resave current" begin
            filename = tempname() * ".rst"
            write(filename, gdalarray)
            gdalarray2 = GDALarray(filename)
            write(gdalarray2)
            @test convert(GeoArray, GDALarray(filename)) == convert(GeoArray, gdalarray2)
        end

        @testset "to grd" begin
            write("testgrd.gri", gdalarray)
            grdarray = GRDarray("testgrd.gri")
            @test crs(grdarray) == convert(ProjString, crs(gdalarray))
            @test bounds(grdarray) == (bounds(gdalarray))
            @test val(dims(grdarray, Y)) == reverse(val(dims(gdalarray, Y)))
            @test val(dims(grdarray, X)) ≈ val(dims(gdalarray, X))
            @test all(GeoArray(grdarray) .== GeoArray(gdalarray))
        end

        @testset "from GeoArray" begin
            filename = tempname() * ".tiff"
            ga = GeoArray(rand(100, 200), (X, Y))
            reorder(ga, (X(ForwardIndex()), Y(ReverseIndex())))
            write(filename, ga)
            saved = GDALarray(filename)
            @test parent(saved[Band(1)]) == parent(ga)
            @test saved[1, 1, 1] == saved[At(1), At(1), At(1)] 
            @test saved[100, 200, 1] == saved[At(100), At(200), At(1)] 
            filename2 = tempname() * ".tif"
            ga2 = GeoArray(rand(100, 200), (X(101:200; mode=Sampled()), Y(1:200; mode=Sampled())))
            write(filename2, ga2)
            @test parent(reorder(GDALarray(filename2)[Band(1)], ForwardArray)) == ga2
        end
       
        # This needs netcdf bounds variables to work
        @testset "to netcdf" begin
            filename2 = tempname() * ".nc"
            write(filename2, gdalarray[Band(1)])
            saved = GeoArray(NCDarray(filename2; crs=crs(gdalarray)))
            @test size(saved) == size(gdalarray[Band(1)])
            @test saved ≈ reverse(gdalarray[Band(1)]; dims=Lat)
            clat, clon = DimensionalData.shiftlocus.(Ref(Center()), dims(gdalarray, (Lat, Lon)))
            @test mappedindex(clat) ≈ reverse(mappedindex(saved, Lat))
            @test mappedindex(clon) ≈ mappedindex(saved, Lon)
            @test all(mappedbounds(saved, Lon) .≈ mappedbounds(clon))
            @test all(mappedbounds(saved, Lat) .≈ mappedbounds(clat))
            # @test all(
            #           projectedbounds(saved, Lon)
            #           locus(saved, Lon)
            #           first(val(span(saved, Lon)))
            #           first(index(saved, Lon))
            #           first(projectedindex(saved, Lon))
            #           last(projectedindex(saved, Lon))
            #           .≈ 
            #           projectedbounds(gdalarray, Lon)
            #           locus(gdalarray, Lon)
            #           first(projectedindex(gdalarray, Lon))
            #           last(projectedindex(gdalarray, Lon))
            #          )
            # @test projectedindex(gdalarray, Lat) .- reverse(projectedindex(DimensionalData.shiftloci(Start(), dims(saved, Lat))))
            @test projectedindex(clon) ≈ projectedindex(saved, Lon)
            # For some reason this crs conversion is less accrurate than the others
            @test all(map((a, b) -> isapprox(a, b; rtol=1e-6), projectedbounds(saved, Lat),  projectedbounds(gdalarray, Lat)))
        end

    end

    @testset "show" begin
        sh = sprint(show, gdalarray)
        # Test but don't lock this down too much
        @test occursin("GDALarray", sh)
        @test occursin("Y", sh)
        @test occursin("X", sh)
        @test occursin("Band", sh)
    end

    @testset "plot and show" begin # TODO write some tests for this
        gdalarray |> show
        gdalarray |> plot
        gdalarray[Y(1)] |> plot
    end

end

@testset "GDAL stack" begin
    gdalstack = stack((a=path, b=path));

    @testset "child array properties" begin
        @test size(gdalstack[:a]) == (514, 515, 1)
        @test gdalstack[:a] isa GeoArray{UInt8,3}
    end

    @testset "indexing" begin
        @test gdalstack[:a][Y(2:3), X(1), Band(1)] == [0x00, 0x6b]
        @test gdalstack[:a][Y(1), X(1), Band(1)] == 0x00
        @test gdalstack[:b, Band(1)] == gdalstack[:b][Band(1)]
        @test typeof(gdalstack[:b, Band(1)]) == typeof(gdalstack[:b][Band(1)])
    end

    @testset "window" begin
        windowedstack = GDALstack((a=path, b=path); window=(Y(1:5), X(1:5), Band(1)))
        @test window(windowedstack) == (Y(1:5), X(1:5), Band(1))
        windowedarray = GeoArray(windowedstack[:a])
        @test windowedarray isa GeoArray{UInt8,2}
        @test length.(dims(windowedarray)) == (5, 5)
        @test size(windowedarray) == (5, 5)
        @test windowedarray[1:3, 2:2] == reshape([0x00, 0x00, 0x00], 3, 1)
        @test windowedarray[1:3, 2] == [0x00, 0x00, 0x00]
        @test windowedarray[1, 2] == 0x00
        windowedstack = GDALstack((a=path, b=path); window=(Y(1:5), X(1:5), Band(1:1)))
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
            geoA = zero(GeoArray(gdalstack[:a]))
            copy!(geoA, gdalstack, :a)
            # First wrap with GeoArray() here or == loads from disk for each cell.
            # we need a general way of avoiding this in all disk-based sources
            @test all(geoA .== GeoArray(gdalstack[:a]))
        end
    end

    @testset "save" begin
        geoA = GeoArray(gdalstack[:a])
        filename = tempname() * ".tif"
        write(filename, gdalstack)
        base, ext = splitext(filename)
        filename_b = string(base, "_b", ext)
        saved = GeoArray(GDALarray(filename_b))
        @test all(saved .== geoA)
    end

end

@testset "GDAL series" begin
    ser = series([path, path], (Ti(),); childkwargs=(mappedcrs=EPSG(4326), name=:test))
    @test GeoArray(ser[Ti(1)]) == GeoArray(GDALarray(path; mappedcrs=EPSG(4326), name=:test))

    gdalstack = GDALstack((a=path, b=path); childtype=GDALarray, childkwargs=(mappedcrs=EPSG(4326),))
    ser = GeoSeries([gdalstack, gdalstack], (Ti,))
    @test ser[1].childkwargs == gdalstack.childkwargs
    # Rebuild the ser by wrapping the GDALarray data in Array.
    # `modify` forces `rebuild` on all containers as in-Memory variants
    modified_ser = modify(Array, ser)
    @test typeof(modified_ser) <: GeoSeries{<:GeoStack{<:NamedTuple{(:a,:b),<:Tuple{<:GeoArray{UInt8,3,<:Tuple,<:Tuple,<:Array{UInt8,3}},Vararg}}}}
end

nothing
