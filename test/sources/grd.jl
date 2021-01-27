using GeoData, Test, Statistics, Dates, Plots
import NCDatasets, ArchGDAL
using GeoData: name, mode, window, DiskStack
testpath = joinpath(dirname(pathof(GeoData)), "../test/")
include(joinpath(testpath, "test_utils.jl"))

maybedownload("https://raw.githubusercontent.com/rspatial/raster/master/inst/external/rlogo.grd", "rlogo.grd")
maybedownload("https://github.com/rspatial/raster/raw/master/inst/external/rlogo.gri", "rlogo.gri")
path = joinpath(testpath, "data/rlogo")
@test isfile(path * ".grd")
@test isfile(path * ".gri")

@testset "Grd array" begin
    grdarray = GRDarray(path);

    @testset "open" begin
        @test all(open(A -> A[Lat=1], grdarray) .=== grdarray[:, 1, :])
    end

    @testset "array properties" begin
        @test grdarray isa GRDarray{Float32,3}
    end

    @testset "dimensions" begin
        @test length(val(dims(dims(grdarray), Lon))) == 101
        @test ndims(grdarray) == 3
        @test dims(grdarray) isa Tuple{<:Lon,<:Lat,<:Band}
        @test refdims(grdarray) == ()
        @test bounds(grdarray) == ((0.0, 101.0), (0.0, 77.0), (1, 3))
    end

    @testset "other fields" begin
        @test missingval(grdarray) == -3.4f38
        @test metadata(grdarray) isa GRDarrayMetadata
        @test name(grdarray) == Symbol("red:green:blue")
        @test label(grdarray) == "red:green:blue"
        @test units(grdarray) == nothing
        customgrdarray = GRDarray(path; name=:test, mappedcrs=EPSG(4326));
        @test name(customgrdarray) == :test
        @test label(customgrdarray) == "test"
        @test mappedcrs(dims(customgrdarray, Lat)) == EPSG(4326)
        @test mappedcrs(dims(customgrdarray, Lon)) == EPSG(4326)
        @test mappedcrs(customgrdarray) == EPSG(4326)
        proj = ProjString("+proj=merc +datum=WGS84")
        @test crs(dims(customgrdarray, Lat)) == proj
        @test crs(dims(customgrdarray, Lon)) == proj
        @test crs(customgrdarray) == proj
    end

    @testset "getindex" begin
        @test grdarray[Band(1)] isa GeoArray{Float32,2}
        @test grdarray[Lat(1), Band(1)] isa GeoArray{Float32,1}
        @test grdarray[Lon(1), Band(1)] isa GeoArray{Float32,1}
        @test grdarray[Lon(50), Lat(30), Band(1)] == 115.0f0
        @test grdarray[1, 1, 1] == 255.0f0
        @test grdarray[Lat(At(20)), Lon(At(20)), Band(3)] == 255.0f0
        @test grdarray[Lat(Contains(60)), Lon(Contains(20)), Band(1)] == 255.0f0
    end

    # @testset "setindex" begin
    #     A = grdarray[:, :, :]
    #     temp = grdarray[1, 1, 1]
    #     println(temp)
    #     @test temp != 100.0f0
    #     grdarray[1, 1, 1] = 100.0f0
    #     grdarray[:, :, :] = 100.0f0
    #     @test grdarray[1, 1, 1] == 100.0f0
    #     grdarray[1, 1, 1] = temp
    #     @test grdarray[1, 1, 1] == temp
    #     println("sum: ", sum(A .- grdarray[:, :, :]))
    #     temp = grdarray[Lon(20), Lat(10), Band(3)]
    #     println(temp)
    #     @test temp != 200.0f0
    #     grdarray[Lon(20), Lat(10), Band(3)] = 200.0f0
    #     @test grdarray[20, 10, 3] == 200.0f0
    #     grdarray[Lon(20), Lat(10), Band(3)] = temp
    # end

    @testset "selectors" begin
        geoarray = grdarray[Lat(Contains(3)), Lon(:), Band(1)]
        @test geoarray isa GeoArray{Float32,1}
        @test grdarray[Lon(Contains(20)), Lat(Contains(10)), Band(1)] isa Float32
    end

    @testset "conversion to GeoArray" begin
        geoarray = grdarray[Lon(1:50), Lat(1:1), Band(1)]
        @test size(geoarray) == (50, 1)
        @test eltype(geoarray) <: Float32
        @time geoarray isa GeoArray{Float32,1}
        @test dims(geoarray) isa Tuple{<:Lon,Lat}
        @test refdims(geoarray) isa Tuple{<:Band}
        @test metadata(geoarray) == metadata(grdarray)
        @test missingval(geoarray) == -3.4f38
        @test name(geoarray) == Symbol("red:green:blue")
    end

    @testset "save" begin
        @testset "2d" begin
            filename2 = tempname()
            write(filename2, GRDarray, grdarray[Band(1)])
            saved = GeoArray(GRDarray(filename2))
            # 1 band is added again on save
            @test size(saved) == size(grdarray[Band(1:1)])
            @test data(saved) == data(grdarray[Band(1:1)])
        end

        @testset "3d with subset" begin
            geoarray = GeoArray(grdarray)[1:100, 1:50, 1:2]
            filename = tempname()
            write(filename, GRDarray, geoarray)
            saved = GeoArray(GRDarray(filename))
            @test size(saved) == size(geoarray)
            @test refdims(saved) == ()
            @test bounds(saved) == bounds(geoarray)
            @test size(saved) == size(geoarray)
            @test missingval(saved) === missingval(geoarray)
            @test metadata(saved) != metadata(geoarray)
            @test metadata(saved)["creator"] == "GeoData.jl"
            @test all(metadata.(dims(saved)) .== metadata.(dims(geoarray)))
            @test name(saved) == name(geoarray)
            @test all(mode.(dims(saved)) .== mode.(dims(geoarray)))
            @test dims(saved) isa typeof(dims(geoarray))
            @test all(val.(dims(saved)) .== val.(dims(geoarray)))
            @test all(mode.(dims(saved)) .== mode.(dims(geoarray)))
            @test all(metadata.(dims(saved)) .== metadata.(dims(geoarray)))
            @test dims(saved) == dims(geoarray)
            @test all(data(saved) .=== data(geoarray))
            @test saved isa typeof(geoarray)
            @test data(saved) == data(geoarray)
        end

        @testset "to netcdf" begin
            filename2 = tempname()
            write(filename2, NCDarray, grdarray[Band(1)])
            saved = GeoArray(NCDarray(filename2; crs=crs(grdarray)))
            @test size(saved) == size(grdarray[Band(1)])
            @test replace_missing(saved, missingval(grdarray)) ≈ reverse(grdarray[Band(1)]; dims=Lat)
            @test replace_missing(saved, missingval(grdarray)) ≈ reverse(grdarray[Band(1)]; dims=Lat)
            @test index(saved, Lon) ≈ index(grdarray, Lon) .+ 0.5
            @test index(saved, Lat) ≈ index(grdarray, Lat) .+ 0.5
            @test bounds(saved, Lat) == bounds(grdarray, Lat)
            @test bounds(saved, Lon) == bounds(grdarray, Lon)
        end

        @testset "to gdal" begin
            # No Band
            gdalfilename = tempname() * ".tif"
            write(gdalfilename, GDALarray, grdarray[Band(1)])
            size(grdarray)
            gdalarray = GDALarray(gdalfilename)
            # @test convert(ProjString, crs(gdalarray)) == convert(ProjString, EPSG(4326))
            @test val(dims(gdalarray, Lon)) ≈ val(dims(grdarray, Lon))
            @test reverse(val(dims(gdalarray, Lat))) ≈ val(dims(grdarray, Lat))
            @test GeoArray(gdalarray) ≈ permutedims(grdarray[Band(1)], [Lon(), Lat()])
            # 3 Bands
            gdalfilename2 = tempname() * ".tif"
            write(gdalfilename2, GDALarray, grdarray)
            gdalarray2 = GDALarray(gdalfilename2)
            @test all(GeoArray(gdalarray2) .== GeoArray(grdarray))
            @test val(dims(gdalarray2, Band)) == 1:3
        end

    end

    @testset "show" begin
        sh = sprint(show, grdarray)
        # Test but don't lock this down too much
        @test occursin("GRDarray", sh)
        @test occursin("Latitude", sh)
        @test occursin("Longitude", sh)
        @test occursin("Band", sh)
    end

    @testset "plot" begin
        grdarray |> plot
        grdarray[Band(1)] |> plot
    end

end

@testset "Grd stack" begin
    grdstack = GRDstack((a=path, b=path))

    @testset "indexing" begin
        @test grdstack[:a][Lat(20), Lon(20), Band(3)] == 70.0f0
        @test grdstack[:a][Lat([2,3]), Lon(40), Band(2)] == [240.0f0, 246.0f0]
    end

    @testset "child array properties" begin
        @test size(grdstack[:a]) == size(GeoArray(grdstack[:a])) == (101, 77, 3)
        @test grdstack[:a] isa GeoArray{Float32,3}
    end

    @testset "window" begin
        windowedstack = GRDstack((a=path, b=path); window=(Lat(1:5), Lon(1:5), Band(1)))
        @test window(windowedstack) == (Lat(1:5), Lon(1:5), Band(1))
        windowedarray = windowedstack[:a]
        @test windowedarray isa GeoArray{Float32,2}
        @test length.(dims(windowedarray)) == (5, 5)
        @test size(windowedarray) == (5, 5)
        @test windowedarray[1:3, 2:2] == reshape([255.0f0, 255.0f0, 255.0f0], 3, 1)
        @test windowedarray[1:3, 2] == [255.0f0, 255.0f0, 255.0f0]
        @test windowedarray[1, 2] == 255.0f0
        windowedstack = GRDstack((a=path, b=path); window=(Lat(1:5), Lon(1:5), Band(1:1)))
        windowedarray = windowedstack[:b]
        @test windowedarray[1:3, 2:2, 1:1] == reshape([255.0f0, 255.0f0, 255.0f0], 3, 1, 1)
        @test windowedarray[1:3, 2:2, 1] == reshape([255.0f0, 255.0f0, 255.0f0], 3, 1)
        @test windowedarray[1:3, 2, 1] == [255.0f0, 255.0f0, 255.0f0]
        @test windowedarray[1, 2, 1] == 255.0f0
        windowedstack = GRDstack((a=path, b=path); window=(Band(1),));
        windowedarray = windowedstack[:b]
        @test windowedarray[1:3, 2:2] == reshape([255.0f0, 255.0f0, 255.0f0], 3, 1)
        @test windowedarray[1:3, 2] == [255.0f0, 255.0f0, 255.0f0]
        @test windowedarray[30, 30] == 185.0f0
        windowedarray |> plot
    end

    # Stack Constructors
    @testset "conversion to GeoStack" begin
        stack = GeoStack(grdstack)
        @test Symbol.(Tuple(keys(grdstack))) == keys(stack)
        smallstack = GeoStack(grdstack; keys=(:a,))
        @test keys(smallstack) == (:a,)
    end

    if VERSION > v"1.1-"
        @testset "copy" begin
            geoarray = zero(GeoArray(grdstack[:a]))
            copy!(geoarray, grdstack, :a)
            # First wrap with GeoArray() here or == loads from disk for each cell.
            # we need a general way of avoiding this in all disk-based sources
            @test geoarray == GeoArray(grdstack[:a])
        end
    end

    @testset "save" begin
        geoarray = GeoArray(grdstack[:b])
        filename = tempname()
        write(filename, GRDarray, grdstack)
        base, ext = splitext(filename)
        filename_b = string(base, "_b", ext)
        saved = GeoArray(GRDarray(filename_b))
        @test typeof(saved) == typeof(geoarray)
        @test data(saved) == data(geoarray)
    end

end

@testset "Grd series" begin
    series = GeoSeries([path, path], (Ti,); childtype=GRDarray, childkwargs=(mappedcrs=EPSG(4326), name=:test))
    @test GeoArray(series[Ti(1)]) ==
        GeoArray(GRDarray(path; mappedcrs=EPSG(4326), name=:test))
    stacks = [DiskStack((a=path, b=path); childtype=GRDarray, childkwargs=(mappedcrs=EPSG(4326), name=:test))]
    series = GeoSeries(stacks, (Ti,))
    @test series[Ti(1)][:a] ==
        GeoArray(GRDarray(path; mappedcrs=EPSG(4326), name=:test))
    modified_series = modify(Array, series)
    @test typeof(modified_series) <: GeoSeries{<:GeoStack{<:NamedTuple{(:a,:b),<:Tuple{<:GeoArray{Float32,3,<:Tuple,<:Tuple,<:Array{Float32,3}},Vararg}}}}
end


nothing
