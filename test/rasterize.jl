using Rasters, Test, ArchGDAL, ArchGDAL.GDAL, Dates, Statistics, DataFrames, Extents, Shapefile, GeometryBasics
import GeoInterface
using Rasters.LookupArrays, Rasters.Dimensions 
using Rasters: bounds

include(joinpath(dirname(pathof(Rasters)), "../test/test_utils.jl"))

A = [missing 7.0f0; 2.0f0 missing]
B = [1.0 0.4; 2.0 missing]
ga = Raster(A, (X, Y); missingval=missing) 
st = RasterStack((a=A, b=B), (X, Y); missingval=(a=missing,b=missing))

pointvec = [(-20.0, 30.0),
            (-20.0, 10.0),
            (0.0, 10.0),
            (0.0, 30.0),
            (-20.0, 30.0)]
vals = [1, 2, 3, 4, 5]
polygon = ArchGDAL.createpolygon(pointvec)
multi_polygon = ArchGDAL.createmultipolygon([[pointvec]])
multi_polygon = ArchGDAL.createmultipolygon([[pointvec]])
multi_point = ArchGDAL.createmultipoint(pointvec)
linestring = ArchGDAL.createlinestring(pointvec)
multi_linestring = ArchGDAL.createmultilinestring([pointvec])
linearring = ArchGDAL.createlinearring(pointvec)
pointfc = map(GeoInterface.getpoint(polygon), vals) do geom, v
    (geometry=geom, val1=v, val2=2.0f0v)
end
table = (X=first.(pointvec), Y=last.(pointvec), othercol=zero.(last.(pointvec)))

test_shape_dir = realpath(joinpath(dirname(pathof(Shapefile)), "..", "test", "shapelib_testcases"))
shp_paths = filter(x -> occursin("shp", x), readdir(test_shape_dir; join=true))
shppath = shp_paths[1]
shphandle = Shapefile.Handle(shppath)

A1 = Raster(zeros(X(-20:5; sampling=Intervals()), Y(0:30; sampling=Intervals())))
A2 = Raster(zeros(Y(0:30; sampling=Intervals()), X(-20:5; sampling=Intervals())))
st = RasterStack((A1, copy(A1)))

@testset "all geoms work as :point" begin
    A = A2
    geom = polygon
    
    for A in (A1, A2), geom in (pointvec, pointfc, multi_point, linestring, multi_linestring, linearring, polygon, multi_polygon, table)
        fill!(A, 0)
        rasterize!(last, A, geom; shape=:point, fill=1);
        @test sum(A) == 4
        fill!(A, 0)
        if !Tables.istable(geom)
            rasterize!(sum, A, [geom, geom]; shape=:point, fill=1)
            @test sum(A) == 10
            A .= 0
        end
    end
    for A in (A1, A2), geom in (table, pointvec, pointfc, multi_point, linestring, multi_linestring, linearring, polygon, multi_polygon)
        st.layer1 .= st.layer2 .= 0
        rasterize!(last, st, geom; shape=:point, fill=(layer1=2, layer2=3))
        @test sum(st[:layer1]) == 8
        @test sum(st[:layer2]) == 12
        @test parent(st[:layer1]) isa Array{Float64,2}
        st[:layer1] .= 0
        st[:layer2] .= 0
        @test_nowarn rasterize!(sum, A, geom; shape=:point, fill=1)
    end
end

@testset "all line and polygon geoms work as :line" begin
    A = A1
    geom = linestring
    for A in (A1, A2), geom in (linestring, multi_linestring, linearring, polygon, multi_polygon)
        A .= 0
        rasterize!(A, geom; shape=:line, fill=1)
        @test sum(A) == 20 + 20 + 20 + 20
    end
    @testset ":line is detected for line geometries" begin
        for A in (A1, A2), geom in (linestring, multi_linestring)
            A .= 0
            rasterize!(A, geom; fill=1)
            @test sum(A) == 20 + 20 + 20 + 20
        end
    end
end

@testset "polygon geoms work as :polygon" begin
    A = A2
    poly = multi_polygon
    for A in (A1, A2), poly in (polygon, multi_polygon)
        ra = rasterize(poly; to=A, missingval=0, shape=:polygon, fill=1, boundary=:center)
        ra_res = rasterize(poly; res=map(step, span(A)), missingval=0, shape=:polygon, fill=1, boundary=:center)
        @test parent(ra) isa Array{Int,2}
        @test sum(ra) == sum(ra_res) === 20 * 20
        ra = rasterize(poly; to=A, shape=:polygon, fill=1, boundary=:touches)
        @test parent(ra) isa Array{Union{Missing,Int},2}
        @test sum(skipmissing(ra)) === 21 * 21
        rasterize!(A, poly; shape=:polygon, fill=1, boundary=:inside)
        @test sum(A) === 19.0 * 19.0
        A .= 0

        @testset "polygon is detected for polygon geometries" begin
            A = Raster(zeros(X(-20:5; sampling=Intervals()), Y(0:30; sampling=Intervals())))
            R = rasterize(poly; to=A, fill=1)
            @test parent(R) isa Array{}
            @test sum(skipmissing(R)) == 20 * 20
            @test_throws ArgumentError rasterize!(A, poly; shape=:notashape, fill=1)
            @test_throws ArgumentError rasterize!(A, poly; shape=:polygon, fill=1, boundary=:notaboundary)
        end

        st1 = rasterize(poly; fill=(layer1=1, layer2=2), to=st)
        @test sum(skipmissing(st1[:layer1])) == 400 # The last value overwrites the first
        @test sum(skipmissing(st1[:layer2])) == 800
        # Missing size / res
        @test_throws ArgumentError rasterize(poly; fill=1)
        # Both size + res
        @test_throws ArgumentError rasterize(poly; res=0.1, size=200, fill=1)
        @test_throws ArgumentError rasterize(poly; res=(0.1, 0.2), size=200, fill=1)
        @test_throws ArgumentError rasterize(poly; res=0.1, size=(200, 200), fill=1)
        @test_throws ArgumentError rasterize(poly; res=(0.1, 0.2), size=(200, 200), fill=1)
    end
end

@testset "from geometries, tables and features of points" begin
    A = A1
    data = DataFrame(pointfc)
    data = multi_point
    data = pointfc

    for data in (pointfc, DataFrame(pointfc), multi_point, pointvec, reverse(pointvec))
        @test sum(skipmissing(rasterize(data; to=A, fill=1))) == 4
        @testset "to and fill Keywords are required" begin
            @test_throws ArgumentError R = rasterize(data; fill=1) 
            @test_throws UndefKeywordError R = rasterize(data; to=A) 
        end
        @testset "NamedTuple of value fill makes a stack" begin
            rst = rasterize(data; to=A, fill=(fill1=3, fill2=6.0f0))
            @test keys(rst) == (:fill1, :fill2)
            @test dims(rst) == dims(A)
            @test map(sum ∘ skipmissing, rst) === (fill1=12, fill2=24.0f0)
        end
        @testset "Single value fill makes an array (ignoring table vals)" begin
            ra = rasterize(data; to=A, fill=0x03, missingval=0x00)
            @test eltype(ra) == UInt8
            @test sum(ra) == 12
        end
    end

    @testset "a single feature" begin
        feature = pointfc[4]
        GeoInterface.isfeature(feature)
        @testset "NTuple of Symbol fill makes an stack" begin
            rst = rasterize(feature; to=A, fill=(:val1, :val2))
            @test parent(rst.val1) isa Array{Union{Missing,Int},2}
            @test parent(rst.val2) isa Array{Union{Missing,Float32},2}
            @test keys(rst) == (:val1, :val2)
            @test dims(rst) == dims(A)
            @test map(sum ∘ skipmissing, rst) === (val1=4, val2=8.0f0)
        end
        @testset "Symbol fill makes an array" begin
            ra = rasterize(feature; to=A, fill=:val1)
            @test ra isa Raster
            @test parent(ra) isa Array{Union{Missing,Int},2}
            @test name(ra) == :val1
        end
    end

    @testset "feature collection, table from fill of Symbol keys" begin
        data = pointfc
        Tables.istable(data)
        for data in (pointfc, DataFrame(pointfc))
            @testset "NTuple of Symbol fill makes an stack" begin
                rst = rasterize(data; to=A, fill=(:val1, :val2))
                @test parent(rst.val1) isa Array{Union{Missing,Int},2}
                @test parent(rst.val2) isa Array{Union{Missing,Float32},2}
                @test keys(rst) == (:val1, :val2)
                @test map(sum ∘ skipmissing, rst) === (val1=14, val2=28.0f0)
                @test_throws ArgumentError rasterize(data; to=A, fill=(:val1, :not_a_column))
            end
            @testset "Symbol fill makes an array" begin
                ra = rasterize(data; to=A, fill=:val1)
                @test parent(ra) isa Array{Union{Missing,Int},2}
                @test ra isa Raster
                @test name(ra) == :val1
                @test_throws ArgumentError rasterize(data; to=A, fill=:not_a_column)
            end
        end
    end
end

function gdal_read_rasterize(fn, options...)
    source_ds = ArchGDAL.unsafe_read(fn)
    dest_ds = ArchGDAL.Dataset(GDAL.gdalrasterize(
        "gadm36_ARG_1.mem",
        Ptr{GDAL.GDALDatasetH}(C_NULL),
        source_ds.ptr,
        GDAL.gdalrasterizeoptionsnew([
            "-of", "MEM",
            "-ts", "250", "250",  # target size:
            "-ot", "Byte",
            "-burn", "1.0",
            options...
        ], C_NULL),
        C_NULL
    ))
    raster = ArchGDAL.read(dest_ds)
    ArchGDAL.destroy(dest_ds)
    ArchGDAL.destroy(source_ds)
    return raster
end

@testset "center in polygon rasterization" begin
    @time gdal_raster = gdal_read_rasterize(shppath);
    @time rasters_raster = rasterize(shphandle.shapes; 
         size=(250, 250), fill=UInt8(1), missingval=UInt8(0),
    );
    # using Plots
    # heatmap(parent(parent(rasters_raster)))
    # heatmap(reverse(gdal_raster[:, :, 1]; dims=2))
    # Same results as GDAL
    @test sum(gdal_raster) == sum(rasters_raster)
    @test reverse(gdal_raster[:, :, 1]; dims=2) == rasters_raster
end

@testset "line touches rasterization" begin
    gdal_touches_raster = gdal_read_rasterize(shppath, "-at")
    rasters_touches_raster = rasterize(shphandle.shapes; 
        size=(250, 250), fill=UInt64(1), missingval=UInt64(0), boundary=:touches
    )
    # Not quite the same answer as GDAL
    @test_broken sum(gdal_touches_raster) == sum(rasters_touches_raster)
    @test_broken reverse(gdal_touches_raster[:, :, 1], dims=2) == rasters_touches_raster
    @test Int(sum(gdal_touches_raster)) == Int(sum(rasters_touches_raster)) - 2
    # Two pixels differ in the angled line, top right
    # using Plots
    # Plots.heatmap(reverse(gdal_touches_raster[:, :, 1], dims=2))
    # Plots.heatmap(parent(parent(rasters_touches_raster)))

    line = LineString([Point(1.00, 4.50), Point(4.75, 0.75)])
    r1 = Raster(zeros(Bool, X(0.5:1.0:6.5; sampling=Intervals()), Y(0.5:1.0:6.5; sampling=Intervals())));
    r1r = reverse(r1; dims=X) 
    Rasters.rasterize!(r1, line; fill=true)
    Rasters.rasterize!(r1r, line; fill=true)
    r2 = Raster(zeros(Bool, X(0.5:1.00001:6.5; sampling=Intervals()), Y(0.5:1.00001:6.5; sampling=Intervals())))
    r2r = reverse(r2; dims=X) 
    Rasters.rasterize!(r2, line; fill=true)
    Rasters.rasterize!(r2r, line; fill=true)
    r3 = Raster(zeros(Bool, X(0.5:0.99999:6.5; sampling=Intervals()), Y(0.5:0.99999:6.5; sampling=Intervals())))
    r3r = reverse(r3; dims=X) 
    Rasters.rasterize!(r3, line; fill=true)
    Rasters.rasterize!(r3r, line; fill=true)

    fill_itr = (x=[1,2,3], y=[4,5,6])
    fill = [vals for (i, vals...) in zip(1:3, values(fill_itr)...)]

    @test reverse(r1r; dims=X) == r1
    @test reverse(r2r; dims=X) == r2
    @test reverse(r3r; dims=X) == r3

    @test sum(r1) == 8
    @test sum(r2) == 9 # The first pixel is larger so the line touches it
    @test sum(r3) == 8

    @test r1 == [
     0 0 0 0 0 0 0
     0 0 0 1 1 0 0
     0 0 1 1 0 0 0
     0 1 1 0 0 0 0
     1 1 0 0 0 0 0
     0 0 0 0 0 0 0
     0 0 0 0 0 0 0]

    @test parent(r2) == [
     0 0 0 0 1 0
     0 0 0 1 1 0
     0 0 1 1 0 0
     0 1 1 0 0 0
     1 1 0 0 0 0
     0 0 0 0 0 0]

    @test r3 == [
     0 0 0 0 0 0 0
     0 0 0 1 1 0 0
     0 0 1 1 0 0 0
     0 1 1 0 0 0 0
     1 1 0 0 0 0 0
     0 0 0 0 0 0 0
     0 0 0 0 0 0 0]
end
# GDAL doesnt do inside / not touching rasterization, so we have no test against GDAL
@testset "line inside rasterization" begin
    @time gdal_raster = gdal_read_rasterize(shppath, "-at")
    rasters_inside_raster = rasterize(shphandle.shapes; 
        size=(250, 250), fill=UInt8(1), missingval=UInt8(0), boundary=:inside
    )
    # using Plots
    # heatmap(parent(parent(rasters_inside_raster)))
end

@testset "reducing rasterization" begin
    pointvec1 = [(-20.0, 30.0),
                (-20.0, 10.0),
                (0.0, 10.0),
                (0.0, 30.0),
                (-20.0, 30.0)]
    pointvec2 = map(p -> (p[1] + 10, p[2] + 10), pointvec1)
    pointvec3 = map(p -> (p[1] + 20, p[2] + 20), pointvec1)
    pointvec4 = map(p -> (p[1] + 30, p[2] + 30), pointvec1)
    polygon = ArchGDAL.createpolygon(pointvec)
    polygons = ArchGDAL.createpolygon.([[pointvec1], [pointvec2], [pointvec3], [pointvec4]])
    # With fill of 1 these are all the same thing
    for f in (last, first, mean, median, maximum, minimum)
        r = rasterize(f, polygons; res=5, fill=1, boundary=:center, crs=EPSG(4326))
        @test parent(r) isa Array{<:Union{Missing,<:Real},2}
        @test sum(skipmissing(r)) == 12 + 12 + 12 + 16
        @test crs(r) == EPSG(4326)
    end
    for f in (last, maximum)
        r = rasterize(last, polygons; res=5, fill=1:4, boundary=:center)
        @test parent(r) isa Array{Union{Missing,Int},2}
        @test sum(skipmissing(r)) == 12 * 1 + 12 * 2 + 12 * 3 + 16 * 4
    end
    for f in (first, minimum)
        r = rasterize(f, polygons; res=5, fill=1:4, boundary=:center)
        @test parent(r) isa Array{Union{Missing,Int},2}
        @test sum(skipmissing(r)) == 16 * 1 + 12 * 2 + 12 * 3 + 12 * 4
    end
    for f in (mean, median)
        r = rasterize(f, polygons; res=5, fill=1:4, boundary=:center)
        @test parent(r) isa Array{Union{Missing,Float64},2}
        @test sum(skipmissing(r)) == 
            (12 * 1 + 8 * 2 + 8 * 3 + 12 * 4) + (4 * 1.5 + 4 * 2.5 + 4 * 3.5)
    end
    prod_r = rasterize(prod, polygons; res=5, fill=1:4, boundary=:center, filename="test.tif")
    prod_r = rasterize(prod, polygons; res=5, fill=1:4, boundary=:center)
    @test sum(skipmissing(prod_r)) == 
        (12 * 1 + 8 * 2 + 8 * 3 + 12 * 4) + (4 * 1 * 2 + 4 * 2 * 3 + 4 * 3 * 4)

    prod_st = rasterize(prod, polygons; res=5, fill=(a=1:4, b=4:-1:1), missingval=missing, boundary=:center)
    @test all(prod_st.a .=== rot180(prod_st.b))
    @test all(prod_r .=== prod_st.a)
    prod_r_m = rasterize(prod, polygons; res=5, fill=1:4, missingval=-1, boundary=:center)
    prod_st_m = rasterize(prod, polygons; res=5, fill=(a=1:4, b=4.0:-1.0:1.0), missingval=(a=-1, b=-1.0), boundary=:center)
    @test all(prod_st_m.a .=== prod_r_m)
    @test all( prod_st_m.b .=== rot180(Float64.(prod_r_m)))

    r = rasterize(last, polygons; res=5, fill=(a=1, b=2), boundary=:center)
    @test all(r.a .* 2 .=== r.b)
    
    reduced_raster_sum_center = rasterize(sum, polygons; res=5, fill=1, boundary=:center)
    reduced_raster_count_center = rasterize(count, polygons; res=5, fill=1, boundary=:center)
    @test name(reduced_raster_sum_center) == :sum
    @test name(reduced_raster_count_center) == :count
    @test sum(skipmissing(reduced_raster_sum_center)) == 
          sum(skipmissing(reduced_raster_count_center)) == 16 * 4
    reduced_raster_sum_touches = rasterize(sum, polygons; res=5, fill=1, boundary=:touches)
    reduced_raster_count_touches = rasterize(count, polygons; res=5, fill=1, boundary=:touches)
    @test name(reduced_raster_sum_touches) == :sum
    @test name(reduced_raster_count_touches) == :count
    # plot(reduced_raster_count_touches)
    # plot(reduced_raster_sum_touches)
    # This is broken because the raster area isn't big enough
    @test_broken sum(skipmissing(reduced_raster_sum_touches)) == 
          sum(skipmissing(reduced_raster_count_touches)) == 25 * 4
    @test sum(skipmissing(reduced_raster_sum_touches)) == 
          sum(skipmissing(reduced_raster_count_touches)) == 25 * 4 - 9
    # The outlines of these plots should exactly mactch, 
    # with three values of 2 on the diagonal
    # using Plots
    # Plots.plot(reduced_raster; clims=(0, 3))
    # Plots.plot!(polygons; opacity=0.3, fillcolor=:black)
    reduced_center = rasterize(sum, polygons; res=5, fill=1, boundary=:center)
    reduced_touches = rasterize(sum, polygons; res=5, fill=1, boundary=:touches)
    reduced_inside = rasterize(sum, polygons; res=5, fill=1, boundary=:inside)
    # Plots.plot(reduced_inside; clims=(0, 3))
    # Plots.plot(reduced_center; clims=(0, 3))
    # Plots.plot(reduced_touches; clims=(0, 3))
    # Plots.plot!(polygons; opacity=0.3, fillcolor=:black)
    # Its not clear what the results here should be - there 
    # are differences between different implementations.
    # Soon we will define the pixel intervals so we don't need
    # arbitrary choices of which lines are touched for :touches
end

@testset "coverage" begin
    @time covsum = coverage(sum, shphandle.shapes; res=1, scale=10);
    @time covunion = coverage(union, shphandle.shapes; res=1, scale=10);
    @test parent(covsum) isa Array{Float64,2}
    @test parent(covunion) isa Array{Float64,2}
    # using Plots
    # plot(covsum; clims=(0, 2))
    # plot(covunion; clims=(0, 2))
    # plot!(shphandle.shapes; opacity=0.2)
    insidecount = rasterize(count, shphandle.shapes; res=1, scale=10, boundary=:inside);
    touchescount = rasterize(count, shphandle.shapes; res=1, scale=10, boundary=:touches);
    # The main polygon should be identical
    @test all(covsum[X=0..120] .=== covunion[X=0..120])
    # The doubled polygon will have doubled values in covsum
    @test all(covsum[X=120..180] .=== covunion[X=120..190] .* 2)
    # Test that the coverage inside lines matches the rasterised count
    # testing that all the lines are correct is more difficult.
    @test all(mask(covsum; with=insidecount) .=== replace_missing(insidecount, 0.0))
    # And test there is nothing outside of the rasterize touches area
    @test all(mask(covsum; with=touchescount) .=== covsum)
    @test !all(mask(covunion; with=insidecount) .=== covunion)
    # TODO test coverage along all the lines is correct somehow
end


