using Rasters, Test, ArchGDAL, ArchGDAL.GDAL, Dates, Statistics, DataFrames, Extents, Shapefile, GeometryBasics
import GeoInterface as GI
using Rasters.Lookups, Rasters.Dimensions 
using Rasters: bounds

include(joinpath(dirname(pathof(Rasters)), "../test/test_utils.jl"))

A = [missing 7.0f0; 2.0f0 missing]
B = [1.0 0.4; 2.0 missing]
ga = Raster(A, (X, Y); missingval=missing) 
st = RasterStack((a=A, b=B), (X, Y); missingval=(a=missing,b=missing))

pointvec = [
    (-20.0, 30.0),
    (-20.0, 10.0),
    (0.0, 10.0),
    (0.0, 30.0),
    (-20.0, 30.0),
]
vals = [1, 2, 3, 4, 5]
polygon = ArchGDAL.createpolygon(pointvec)
multi_polygon = ArchGDAL.createmultipolygon([[pointvec]])
multi_point = ArchGDAL.createmultipoint(pointvec)
linestring = ArchGDAL.createlinestring(pointvec)
multi_linestring = ArchGDAL.createmultilinestring([pointvec])
linearring = ArchGDAL.createlinearring(pointvec)
line_collection = GI.GeometryCollection([multi_linestring])
poly_collection = GI.GeometryCollection([polygon])
pointtable = map(GI.getpoint(polygon), vals) do geom, v
    (geometry=geom, val1=v, val2=2.0f0v)
end
pointfc = map(GI.getpoint(polygon), vals) do geom, v
    GI.Feature(geom; properties=(val1=v, val2=2.0f0v))
end |> GI.FeatureCollection
pointdf = DataFrame(pointtable)
table = (X=first.(pointvec), Y=last.(pointvec), othercol=zero.(last.(pointvec)))

geoms = (pointvec, pointtable, pointfc, pointdf, multi_point, linestring, multi_linestring, linearring, polygon, multi_polygon)#, table)
collections = (line_collection, poly_collection)

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
    geom = linearring
    geom = pointvec
    geom = line_collection
    geom = poly_collection
    geom = pointfc
    threaded = true
    
    for A in (A1, A2), 
        geom in (geoms..., collections...),
        threaded in (true, false)

        fill!(A, 0)
        rasterize!(sum, A, geom; shape=:point, fill=1, threaded);
        @test sum(A) == 5.0
        @test sum(rasterize(sum, geom; to=A, shape=:point, fill=1, missingval=0, threaded)) == 5.0
        @test sum(rasterize(xs -> sum(xs), geom; to=A, shape=:point, fill=1, missingval=0, threaded)) == 5.0
        @test sum(rasterize(xs -> sum(xs), geom; to=A, shape=:point, fill=1, missingval=0, threaded, threadsafe=true)) == 5.0
        rasterize!(last, A, geom; shape=:point, fill=1, threaded);
        @test sum(A) == 4.0
        @test sum(rasterize(last, geom; to=A, shape=:point, fill=1, missingval=0, threaded)) == 4.0
        fill!(A, 0)
        if !(Tables.istable(geom) || GI.isfeaturecollection(geom))
            rasterize!(count, A, [geom; geom]; shape=:point, threaded)
            @test sum(A) == 10.0
            fill!(A, 0)
        end

        # stack
        rasterize!(sum, st, geom; shape=:point, fill=(layer1=2, layer2=3), threaded)
        st.layer1 .= st.layer2 .= 0
        rasterize!(sum, st, geom; shape=:point, fill=(layer1=2, layer2=3), threaded)
        @test sum(st.layer1) == 10
        @test sum(st.layer2) == 15
        @test parent(st.layer1) isa Array{Float64,2}
        st.layer1 .= 0; st.layer2 .= 0
        rasterize!(last, st, geom; shape=:point, fill=(layer1=2, layer2=3), threaded)
        @test sum(st.layer1) == 8
        @test sum(st.layer2) == 12
        st.layer1 .= 0; st.layer2 .= 0
        rasterize!(first, st, geom; shape=:point, fill=(layer1=2, layer2=3), threaded)
        @test sum(st.layer1) == 8
        @test sum(st.layer2) == 12
        st.layer1 .= 0; st.layer2 .= 0
    end

    for A in (A1, A2), 
        geom in geoms,
        threaded in (true, false)

        st[:layer1] .= 0; st[:layer2] .= 0
        rasterize!(sum, st, geom; shape=:point, fill=(layer1=1:5, layer2=6:10), threaded)
        @test sum(st[:layer1]) == sum(1:5)
        @test sum(st[:layer2]) == sum(6:10)
        @test parent(st[:layer1]) isa Array{Float64,2}
        st[:layer1] .= 0; st[:layer2] .= 0
        rasterize!(last, st, geom; shape=:point, fill=(layer1=1:5, layer2=6:10), threaded)
        @test sum(st[:layer1]) == sum(2:5)
        @test sum(st[:layer2]) == sum(7:10)
        st[:layer1] .= 0; st[:layer2] .= 0
        rasterize!(first, st, geom; shape=:point, fill=(layer1=1:5, layer2=6:10), threaded)
        @test sum(st[:layer1]) == sum(1:4)
        @test sum(st[:layer2]) == sum(6:9)
        st[:layer1] .= 0; st[:layer2] .= 0
        @test_nowarn rasterize!(sum, A, geom; shape=:point, fill=1, threaded)
    end
end

@testset "all line and polygon geoms work as :line" begin
    A = A1
    geom = linestring
    geom = line_collection
    for A in (A1, A2), geom in (linestring, multi_linestring, linearring, polygon, multi_polygon, line_collection, poly_collection),
        threaded in (true, false)
        A .= 0
        rasterize!(sum, A, geom; shape=:line, fill=1, threaded)
        @test sum(A) == 20 + 20 + 20 + 20
        @test sum(rasterize(sum, geom; to=A, shape=:line, fill=1, missingval=0, threaded)) == 80
        @test sum(rasterize(xs -> sum(xs), geom; to=A, shape=:line, fill=1, missingval=0, threaded, threadsafe=true)) == 80
        @test sum(rasterize(xs -> sum(xs), geom; to=A, shape=:line, fill=1, missingval=0, threaded, threadsafe=false)) == 80
    end
    @testset ":line is detected for line geometries" begin
        for A in (A1, A2), geom in (linestring, multi_linestring), threaded in (true, false)
            A .= 0
            rasterize!(A, geom; fill=1, threaded)
            @test sum(A) == 20 + 20 + 20 + 20
            @test sum(rasterize(geom; to=A, fill=1, missingval=0, threaded)) == 80
        end
    end
end

@testset "polygon geoms work as :polygon" begin
    A = A1
    poly = polygon
    poly = poly_collection
    threaded = false
    for A in (A1, A2), poly in (polygon, multi_polygon, poly_collection), threaded in (true, false)
        A .= 0
        ra = rasterize(last, poly; to=A, missingval=0, shape=:polygon, fill=1, boundary=:center, threaded)
        ra_res = rasterize(last, poly; res=map(step, span(A)), missingval=0, shape=:polygon, fill=1, boundary=:center, threaded)
        @test parent(ra) isa Matrix{Int}
        @test sum(ra) == sum(ra_res) === 20 * 20
        ra = rasterize(last, poly; to=A, shape=:polygon, fill=1, boundary=:touches, threaded)
        @test parent(ra) isa Array{Union{Missing,Int},2}
        @test sum(skipmissing(ra)) === 21 * 21
        rasterize!(last, A, poly; shape=:polygon, fill=1, boundary=:inside, threaded)
        @test sum(A) === 19.0 * 19.0
        A .= 0

        @testset "polygon is detected for polygon geometries" begin
            A = Raster(zeros(X(-20:5; sampling=Intervals()), Y(0:30; sampling=Intervals())))
            R = rasterize(last, poly; to=A, fill=1, threaded)
            @test parent(R) isa Array{}
            @test sum(skipmissing(R)) == 20 * 20
            @test_throws ArgumentError rasterize!(A, poly; shape=:notashape, fill=1, threaded)
            @test_throws ArgumentError rasterize!(A, poly; shape=:polygon, fill=1, boundary=:notaboundary, threaded)
        end

        st1 = rasterize(last, poly; fill=(layer1=1, layer2=2), to=st, threaded)
        @test sum(skipmissing(st1[:layer1])) == 400 # The last value overwrites the first
        @test sum(skipmissing(st1[:layer2])) == 800
        # Missing size / res
        @test_throws ArgumentError rasterize(poly; fill=1, threaded)
        # Both size + res
        @test_throws ArgumentError rasterize(poly; res=0.1, size=200, fill=1, threaded)
        @test_throws ArgumentError rasterize(poly; res=(0.1, 0.2), size=200, fill=1, threaded)
        @test_throws ArgumentError rasterize(poly; res=0.1, size=(200, 200), fill=1, threaded)
        @test_throws ArgumentError rasterize(poly; res=(0.1, 0.2), size=(200, 200), fill=1, threaded)
    end
end

@testset "from geometries, tables and features of points" begin
    A = A1
    data = DataFrame(pointtable)
    data = multi_point
    data = pointfc

   for data in (pointtable, pointfc, DataFrame(pointtable), multi_point, pointvec, reverse(pointvec))
        @test sum(skipmissing(rasterize(sum, data; to=A, fill=1))) == 5
        @testset "to and fill Keywords are required" begin
            @test_throws ArgumentError R = rasterize(data; fill=1) 
            @test_throws UndefKeywordError R = rasterize(data; to=A) 
        end
        @testset "NamedTuple of value fill makes a stack" begin
            rst = rasterize(sum, data; to=A, fill=(fill1=3, fill2=6.0f0))
            @test eltype(rst) == @NamedTuple{fill1::Union{Missing,Int64}, fill2::Union{Missing,Float32}}
            @test keys(rst) == (:fill1, :fill2)
            @test dims(rst) == dims(A)
            @test map(sum ∘ skipmissing, rst) === (fill1=15, fill2=30.0f0)
        end
        @testset "Single value fill makes an array (ignoring table vals)" begin
            ra = rasterize(sum, data; to=A, fill=0x03, missingval=0x00)
            @test ra isa Raster
            @test eltype(ra) == UInt64
            @test sum(ra) === 0x000000000000000f
        end
    end

    @testset "a single feature" begin
        feature = pointtable[4]
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
        data = pointtable
        data = pointfc
        for data in (pointfc, pointtable, DataFrame(pointtable)), threaded in (true, false)
            @testset "NTuple of Symbol fill makes an stack" begin
                rst = rasterize(sum, data; to=A, fill=(:val1, :val2), threaded)
                @test parent(rst.val1) isa Array{Union{Missing,Int},2}
                @test parent(rst.val2) isa Array{Union{Missing,Float32},2}
                @test keys(rst) == (:val1, :val2)
                @test map(sum ∘ skipmissing, rst) === (val1=15, val2=30.0f0)
                @test_throws ArgumentError rasterize(data; to=A, fill=(:val1, :not_a_column), threaded)
            end
            @testset "Symbol fill makes an array" begin
                ra = rasterize(sum, data; to=A, fill=:val1)
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
    @time rasters_raster = rasterize(last, shphandle.shapes; 
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
    rasters_touches_raster = rasterize(last, shphandle.shapes; 
        size=(250, 250), fill=UInt64(1), missingval=UInt64(0), boundary=:touches
    )
    # Not quite the same answer as GDAL
    @test_broken sum(gdal_touches_raster) == sum(rasters_touches_raster)
    @test_broken reverse(gdal_touches_raster[:, :, 1], dims=2) == rasters_touches_raster
    # Test that its knwon to be off by 2:
    @test count(reverse(gdal_touches_raster[:, :, 1], dims=2) .== rasters_touches_raster) == length(rasters_touches_raster) - 2
    # Two pixels differ in the angled line, top right
    # using Plots
    # Plots.heatmap(reverse(gdal_touches_raster[:, :, 1], dims=2))
    # Plots.heatmap(parent(parent(rasters_touches_raster)))

    line = LineString([Point(1.00, 4.50), Point(4.75, 0.75)])
    r1 = Raster(zeros(Bool, X(0.5:1.0:6.5; sampling=Intervals()), Y(0.5:1.0:6.5; sampling=Intervals())));
    r1r = reverse(r1; dims=X) 
    Rasters.rasterize!(last, r1, line; fill=true)
    Rasters.rasterize!(last, r1r, line; fill=true)
    r2 = Raster(zeros(Bool, X(0.5:1.00001:6.5; sampling=Intervals()), Y(0.5:1.00001:6.5; sampling=Intervals())))
    r2r = reverse(r2; dims=X) 
    Rasters.rasterize!(last, r2, line; fill=true)
    Rasters.rasterize!(last, r2r, line; fill=true)
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
# GDAL doesn't do inside / not touching rasterization, so we have no test against GDAL
@testset "line inside rasterization" begin
    @time gdal_raster = gdal_read_rasterize(shppath, "-at")
    rasters_inside_raster = rasterize(last, shphandle.shapes; 
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
    for  threaded in (true, false)
        for f in (last, first, mean, median, maximum, minimum)
            r = rasterize(f, polygons; res=5, fill=1, boundary=:center, threaded, crs=EPSG(4326))
            @test parent(r) isa Array{<:Union{Missing,<:Real},2}
            @test sum(skipmissing(r)) == 12 + 12 + 12 + 16
            @test crs(r) == EPSG(4326)
        end
        for f in (last, maximum)
            r = rasterize(last, polygons; res=5, fill=1:4, boundary=:center, threaded)
            @test parent(r) isa Array{Union{Missing,Int},2}
            @test sum(skipmissing(r)) == 12 * 1 + 12 * 2 + 12 * 3 + 16 * 4
        end
        for f in (first, minimum)
            r = rasterize(f, polygons; res=5, fill=1:4, boundary=:center, threaded)
            @test parent(r) isa Array{Union{Missing,Int},2}
            @test sum(skipmissing(r)) == 16 * 1 + 12 * 2 + 12 * 3 + 12 * 4
        end
        for f in (mean, median)
            r = rasterize(f, polygons; res=5, fill=1:4, boundary=:center)
            @test parent(r) isa Array{Union{Missing,Float64},2}
            @test sum(skipmissing(r)) == 
                (12 * 1 + 8 * 2 + 8 * 3 + 12 * 4) + (4 * 1.5 + 4 * 2.5 + 4 * 3.5)
        end
        filename = tempname() * ".tif"
        prod_r = rasterize(prod, polygons; res=5, fill=1:4, boundary=:center, filename, threaded)
        prod_r = rasterize(prod, polygons; res=5, fill=1:4, boundary=:center, threaded)
        @test sum(skipmissing(prod_r)) == 
            (12 * 1 + 8 * 2 + 8 * 3 + 12 * 4) + (4 * 1 * 2 + 4 * 2 * 3 + 4 * 3 * 4)

        prod_st = rasterize(prod, polygons; res=5, fill=(a=1:4, b=4:-1:1), missingval=missing, boundary=:center, threaded)
        @test all(prod_st.a .=== rot180(parent(prod_st.b)))
        @test all(prod_r .=== prod_st.a)
        prod_r_m = rasterize(prod, polygons; res=5, fill=1:4, missingval=-1, boundary=:center, threaded)
        prod_st_m = rasterize(prod, polygons; res=5, fill=(a=1:4, b=4.0:-1.0:1.0), missingval=(a=-1, b=-1.0), boundary=:center, threaded)
        @test all(prod_st_m.a .=== prod_r_m)
        @test all(prod_st_m.b .=== rot180(parent(Float64.(prod_r_m))))

        r = rasterize(last, polygons; res=5, fill=(a=1, b=2), boundary=:center, threaded)
        @test all(r.a .* 2 .=== r.b)
        
        reduced_raster_sum_center = rasterize(sum, polygons; res=5, fill=1, boundary=:center, threaded)
        reduced_raster_count_center = rasterize(count, polygons; res=5, fill=1, boundary=:center, threaded)
        @test name(reduced_raster_sum_center) == :sum
        @test name(reduced_raster_count_center) == :count
        @test sum(skipmissing(reduced_raster_sum_center)) == 
              sum(skipmissing(reduced_raster_count_center)) == 16 * 4
        reduced_raster_sum_touches = rasterize(sum, polygons; res=5, fill=1, boundary=:touches, threaded)
        reduced_raster_count_touches = rasterize(count, polygons; res=5, fill=1, boundary=:touches, threaded)
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
        reduced_center = rasterize(sum, polygons; res=5, fill=1, boundary=:center, threaded)
        reduced_touches = rasterize(sum, polygons; res=5, fill=1, boundary=:touches, threaded)
        reduced_inside = rasterize(sum, polygons; res=5, fill=1, boundary=:inside, threaded)
        # Plots.plot(reduced_inside; clims=(0, 3))
        # Plots.plot(reduced_center; clims=(0, 3))
        # Plots.plot(reduced_touches; clims=(0, 3))
        # Plots.plot!(polygons; opacity=0.3, fillcolor=:black)
        # Its not clear what the results here should be - there 
        # are differences between different implementations.
        # Soon we will define the pixel intervals so we don't need
        # arbitrary choices of which lines are touched for :touches
    end
end

@testset "Rasterizing feature collection of polygons" begin
    p1 = GI.Polygon([[[-55965.680060140774, -31588.16072168928], [-55956.50771556479, -31478.09258677756], [-31577.548550575284, -6897.015828572996], [-15286.184961223798, -15386.952072224134], [-9074.387601621409, -27468.20712382156], [-8183.4538916097845, -31040.003969070774], [-27011.85123029944, -38229.02388009402], [-54954.72822634951, -32258.9734800704], [-55965.680060140774, -31588.16072168928]]])
    p2 = GI.Polygon([[[-80000.0, -80000.0], [-80000.0, 80000.0], [-60000.0, 80000.0], [-50000.0, 40000.0], [-60000.0, -80000.0], [-80000.0, -80000.0]]])
    fc = GI.FeatureCollection([GI.Feature(p1, properties = (;val=1)), GI.Feature(p2, properties = (;val=2))])

    vecint = rasterize(sum, [p1, p2]; size=(100, 100), fill=[1, 2])
    fcint = rasterize(sum, fc; size=(100, 100), fill=[1, 2])
    fcval = rasterize(sum, fc; size=(100, 100), fill=:val)
    @test all(vecint .=== fcval .=== fcint)
end

@testset "Rasterizing with extra dimensions" begin
    A3 = cat(A1, A1, A1; dims=Band)
    Rasters.rasterize!(last, A3, polygon; fill=true)
    @test A3[Band(1)] == A3[Band(2)] == A3[Band(3)]
    @test sum(A3) == 1200
end

@testset "Rasterize an empty polygon" begin
    empty_polygon = ArchGDAL.createpolygon(Tuple{Float64,Float64}[])
    rast = rasterize(empty_polygon; to=A1, fill=1, missingval=0)
    @test sum(rast) == 0
end

@testset "Cant rasterize to a new empty array" begin
    # Its difficult to make a new empty raster here (we possibly could in future?)
    @test_throws ArgumentError rasterize(polygon; to=A1[1:0, 1:0], fill=1, missingval=0)
    # But we just warn in `rasterize!` because the result is still correct
    @test_warn "Destination is empty" rasterize!(A1[1:0, 1:0], polygon; to=A1[1:0, 1:0], fill=1, missingval=0)
end

@testset "coverage" begin
    @time covsum = coverage(sum, shphandle.shapes; threaded=false, res=1, scale=10)
    @time covunion = coverage(union, shphandle.shapes; threaded=false, res=1, scale=10)
    @test parent(covsum) isa Array{Float64,2}
    @test parent(covunion) isa Array{Float64,2}

    @time covsum = coverage(sum, shphandle.shapes; threaded=true, res=1, scale=10)
    @time covunion = coverage(union, shphandle.shapes; threaded=true, res=1, scale=10)
    @test parent(covsum) isa Array{Float64,2}
    @test parent(covunion) isa Array{Float64,2}
    # using Plots
    # plot(covsum; clims=(0, 2))
    # plot(covunion; clims=(0, 2))
    # plot!(shphandle.shapes; opacity=0.2)
    insidecount = rasterize(count, shphandle.shapes; res=1, boundary=:inside);
    touchescount = rasterize(count, shphandle.shapes; res=1, boundary=:touches);
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
    
    @test_throws ErrorException coverage(union, shphandle.shapes; threaded=false, res=1, scale=10000)
    # Too slow and unreliable to test in CI, but it warns and uses one thread given 32gb of RAM: 
    # coverage(union, shphandle.shapes; threaded=true, res=1, scale=1000)
end
