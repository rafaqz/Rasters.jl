using GeoData, Test, Statistics, Dates
using GeoData: name
include("test_utils.jl")


geturl("https://raw.githubusercontent.com/rspatial/raster/master/inst/external/rlogo.grd", "rlogo.grd")
geturl("https://github.com/rspatial/raster/raw/master/inst/external/rlogo.gri", "rlogo.gri")
path = "data/rlogo"
@test isfile(path * ".grd")
@test isfile(path * ".gri")

@testset "array" begin
    grdarray = GeoData.GrdArray(path);

    @testset "array properties" begin
        @test grdarray isa GrdArray{Float32,3}
    end

    @testset "dimensions" begin
        @test length(val(dims(dims(grdarray), Lon))) == 101
        @test ndims(grdarray) == 3
        @test dims(grdarray) isa Tuple{<:Lon,<:Lat,<:Band}
        @test refdims(grdarray) == ()
        @test bounds(grdarray) == ((0.0, 77.0), (0.0, 101.0), (1, 3))
    end

    @testset "other fields" begin
        @test GeoData.window(grdarray) == ()
        @test missingval(grdarray) == -3.4f38 
        @test metadata(grdarray) isa GrdMetadata
        @test name(grdarray) == "red:green:blue"
    end

    @testset "getindex" begin 
        @test grdarray[Band(1)] isa GeoArray{Float32,2} 
        @test grdarray[Lat(1), Band(1)] isa GeoArray{Float32,1} 
        @test grdarray[Lon(1), Band(1)] isa GeoArray{Float32,1}
        @test grdarray[Lon(1), Lat(1), Band(1)] isa Float32 
        @test grdarray[1, 1, 1] isa Float32
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
        geoarray = grdarray[Lat(Near(3)), Lon(:), Band(1)]
        @test geoarray isa GeoArray{Float32,1}
        # @test bounds(a) == ()
        @test grdarray[Lon(Near(20)), Lat(Near(10)), Band(1)] isa Float32
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
        @test name(geoarray) == "red:green:blue"
    end

    @testset "save" begin
        geoarray = grdarray[Band(1)]
        filename = tempname()
        write(filename, GrdArray, geoarray)
        saved = GeoArray(GrdArray(filename))
        # 1 bands is added again on save
        @test size(saved) != size(geoarray)
        @test refdims(saved) == ()
        saved = saved[Band(1)]
        @test bounds(saved) == bounds(geoarray)
        @test size(saved) == size(geoarray)
        @test missingval(saved) === missingval(geoarray)
        @test metadata(saved) != metadata(geoarray)
        @test metadata(saved)["creator"] == "GeoData.jl"
        @test all(metadata.(dims(saved)) .== metadata.(dims(geoarray)))
        @test GeoData.name(saved) == GeoData.name(geoarray)
        @test all(DimensionalData.grid.(dims(saved[Band(1)])) .== DimensionalData.grid.(dims(geoarray)))
        @test dims(saved) isa typeof(dims(geoarray))
        @test all(val.(dims(saved)) .== val.(dims(geoarray)))
        @test all(metadata.(dims(saved)) .== metadata.(dims(geoarray)))
        @test all(data(saved) .=== data(geoarray))
        @test saved isa typeof(geoarray)
        write(filename, GrdArray, grdarray)
        saved = GeoArray(GrdArray(filename))
        @test size(saved) == size(grdarray)
    end

end

@testset "stack" begin
    grdstack = GeoStack((a=GrdArray(path), b=GrdArray(path)))

    @test grdstack[:a][Lat(1), Lon(1), Band(1)] == 255.0f0
    @test grdstack[:a][Lat([2,3]), Lon(1), Band(1)] == [255.0f0, 255.0f0] 

    # Stack Constructors
    @testset "conversion to GeoStack" begin
        stack = GeoStack(grdstack)
        @test Symbol.(Tuple(keys(grdstack))) == keys(stack)
        smallstack = GeoStack(grdstack; keys=(:a,))
        keys(smallstack) == (:a,)
    end

    if VERSION > v"1.1-"
        @testset "copy" begin
            geoarray = zero(GeoArray(grdstack[:a]))
            copy!(geoarray, grdstack, :a)
            maximum(grdstack[:a])
            maximum(geoarray)
            # First wrap with GeoArray() here or == loads from disk for each cell.
            # we need a general way of avoiding this in all disk-based sources
            @test geoarray == GeoArray(grdstack[:a])
        end
    end

    @testset "save" begin
        geoarray = GeoArray(grdstack[:a])
        filename = tempname()
        write(filename, GrdArray, grdstack)
        base, ext = splitext(filename)
        filename_b = string(base, "_b", ext)
        saved = GeoArray(GrdArray(filename_b))
        @test saved == geoarray
    end

end
