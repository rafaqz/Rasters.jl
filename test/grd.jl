using GeoData, Test, Statistics, Dates
include("utils.jl")

grdpath = geturl("https://raw.githubusercontent.com/rspatial/raster/master/inst/external/rlogo.grd")
gripath = geturl("https://github.com/rspatial/raster/blob/master/inst/external/rlogo.gri?raw=true", "rlogo.gri")
path = "rlogo"
@test isfile(path * ".grd")

@testset "array" begin
    grdarray = GeoData.GrdArray(path);

    @testset "array properties" begin
        @test typeof(grdarray) <: GrdArray{Float32,3}
    end

    @testset "dimensions" begin
        @test length(val(dims(dims(grdarray), Lon))) == 101
        @test ndims(grdarray) == 3
        @test typeof(dims(grdarray)) <: Tuple{<:Lon,<:Lat,<:Band}
        @test refdims(grdarray) == ()
        @test bounds(grdarray) == ((0.0, 77.0), (0.0, 101.0), (1, 3))
    end

    @testset "other fields" begin
        @test GeoData.window(grdarray) == ()
        @test missingval(grdarray) == -3.4f38 
        @test typeof(metadata(grdarray)) <: GrdMetadata
        @test name(grdarray) == "red:green:blue"
    end

    @testset "getindex" begin 
        @test typeof(grdarray[Band(1)]) <: GeoArray{Float32,2} 
        @test typeof(grdarray[Lat(1), Band(1)]) <: GeoArray{Float32,1} 
        @test typeof(grdarray[Lon(1), Band(1)]) <: GeoArray{Float32,1}
        @test typeof(grdarray[Lon(1), Lat(1), Band(1)]) <: Float32 
        @test typeof(grdarray[1, 1, 1]) <: Float32
    end

    @testset "setindex" begin 
        # Make a copy to write to
        path2 = "rlogo2"
        cp(path * ".grd", path2 * ".grd"; force=true)
        cp(path * ".gri", path2 * ".gri"; force=true)
        grdarray2 = GeoData.GrdArray(path2);
        @test grdarray2[1, 1, 1] != 100.0f0
        grdarray2[1, 1, 1] = 100.0f0
        @test grdarray2[1, 1, 1] == 100.0f0
        @test grdarray2[20, 10, 3] != 200.0f0
        grdarray2[Lon(20), Lat(10), Band(3)] = 200.0f0
        @test grdarray2[20, 10, 3] == 200.0f0
    end

    @testset "selectors" begin
        geoarray = grdarray[Lat(Near(3)), Lon(:), Band(1)]
        @test typeof(geoarray) <: GeoArray{Float32,1}
        # @test bounds(a) == ()
        # Doesn't handle returning a single value
        x = typeof(grdarray[Lon(Near(20)), Lat(Near(10)), Band(1)]) <: Float32
    end

    @testset "conversion to GeoArray" begin
        geoarray = grdarray[Lon(1:50), Lat(1:1), Band(1)]
        @test size(geoarray) == (50, 1)
        @test eltype(geoarray) <: Float32
        @time typeof(geoarray) <: GeoArray{Float32,1} 
        @test typeof(dims(geoarray)) <: Tuple{<:Lon,Lat}
        @test typeof(refdims(geoarray)) <: Tuple{<:Band} 
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
        @test val(metadata(saved))["creator"] == "GeoData.jl" 
        @test all(metadata.(dims(saved)) .== metadata.(dims(geoarray)))
        @test GeoData.name(saved) == GeoData.name(geoarray)
        @test all(DimensionalData.grid.(dims(saved[Band(1)])) .== DimensionalData.grid.(dims(geoarray)))
        @test typeof(dims(saved)) == typeof(dims(geoarray))
        @test val(dims(saved)[1]) == val(dims(geoarray)[1])
        @test val(dims(saved)[2]) == val(dims(geoarray)[2])
        @test all(val.(dims(saved)) .== val.(dims(geoarray)))
        @test all(metadata.(dims(saved)) .== metadata.(dims(geoarray)))
        @test all(parent(saved) .=== parent(geoarray))
        @test typeof(saved) == typeof(geoarray)
        geoarray = grdarray
        write(filename, GrdArray, geoarray)
        saved = GeoArray(GrdArray(filename))
        @test size(saved) == size(geoarray)
    end

end

@testset "stack" begin
    grdstack = GeoStack((a=GrdArray(path), b=GrdArray(path)))

    @test grdstack[:a][Lat([2,3]), Lon(1), Band(1)] == [255.0, 255.0] 
    @test grdstack[:a][Lat(1), Lon(1), Band(1)] == 255.0

    # Stack Constructors
    @testset "conversion to GeoStack" begin
        stack = GeoStack(grdstack)
        @test Symbol.(Tuple(keys(grdstack))) == keys(stack)
        smallstack = GeoStack(grdstack; keys=(:a,))
        keys(smallstack) == (:a,)
    end

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
