using Rasters, Test, ArchGDAL, ArchGDAL.GDAL, Dates, Statistics, DataFrames, Extents, Shapefile, GeometryBasics
import GeoInterface
using Rasters.LookupArrays, Rasters.Dimensions 
using Rasters: bounds

include(joinpath(dirname(pathof(Rasters)), "../test/test_utils.jl"))

A = [missing 7.0f0; 2.0f0 missing]
B = [1.0 0.4; 2.0 missing]
ga = Raster(A, (X, Y); missingval=missing) 
ga99 = replace_missing(ga, -9999)
gaNaN = replace_missing(ga, NaN32)
gaMi = replace_missing(ga)
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

test_shape_dir = realpath(joinpath(dirname(pathof(Shapefile)), "..", "test", "shapelib_testcases"))
shp_paths = filter(x -> occursin("shp", x), readdir(test_shape_dir; join=true))
shppath = shp_paths[1]
shphandle = Shapefile.Handle(shppath)

@testset "replace_missing" begin
    @test all(isequal.(ga99, [-9999.0f0 7.0f0; 2.0f0 -9999.0f0]))
    @test missingval(ga99) === -9999.0f0
    @test all(isequal.(gaNaN, [NaN32 7.0f0; 2.0f0 NaN32]))
    @test missingval(gaNaN) === NaN32
    @test all(isequal.(gaMi, ga))
    @test missingval(gaMi) === missing
    # The second layer NaN32 will be promoted to NaN to match the array type
    @test missingval(replace_missing(st, NaN32)) === (a=NaN32, b=NaN)
    @test all(map(values(replace_missing(st, NaN32)), (a=[NaN32 7.0f0; 2.0f0 NaN32], b=[1.0 0.4; 2.0 NaN])) do x, y
        all(x .=== y)
    end)
    dNaN = replace_missing(ga, NaN32; filename="test.tif")
    @test all(isequal.(dNaN, [NaN32 7.0f0; 2.0f0 NaN32]))
    rm("test.tif")
    stNaN = replace_missing(st, NaN32; filename="teststack.tif")
    @test all(map(stNaN[Band(1)], (a=[NaN32 7.0f0; 2.0f0 NaN32], b=[1.0 0.4; 2.0 NaN])) do x, y
        all(x .=== y)
    end)
    rm("teststack_a.tif")
    rm("teststack_b.tif")
end

@testset "boolmask" begin
    @test boolmask(ga) == [false true; true false]
    @test boolmask(ga99) == [false true; true false]
    @test boolmask(gaNaN) == [false true; true false]
    @test dims(boolmask(ga)) == (X(NoLookup(Base.OneTo(2))), Y(NoLookup(Base.OneTo(2))))
    @test boolmask(polygon; res=1.0) == trues(X(Projected(-20:1.0:-1.0; crs=nothing)), Y(Projected(10.0:1.0:29.0; crs=nothing)))
    # With a :geometry axis
    x = boolmask([polygon, polygon]; combine=false, res=1.0)
    @test eltype(x) == Bool
    @test size(x) == (20, 20, 2)
    @test sum(x) == 800
    x = boolmask([polygon, polygon]; combine=true, res=1.0)
    @test size(x) == (20, 20)
    @test sum(x) == 400
end

@testset "missingmask" begin
    @test all(missingmask(ga) .=== [missing true; true missing])
    @test all(missingmask(ga99) .=== [missing true; true missing])
    @test all(missingmask(gaNaN) .=== [missing true; true missing])
    @test dims(missingmask(ga)) == (X(NoLookup(Base.OneTo(2))), Y(NoLookup(Base.OneTo(2))))
    @test missingmask(polygon; res=1.0) == fill!(Raster{Union{Missing,Bool}}(undef, X(Projected(-20:1.0:-1.0; crs=nothing)), Y(Projected(10.0:1.0:29.0; crs=nothing))), true)
    x = missingmask([polygon, polygon]; combine=false, res=1.0)
    @test eltype(x) == Union{Bool,Missing}
    @test size(x) == (20, 20, 2)
    @test sum(x) == 800
    x = missingmask([polygon, polygon]; combine=true, res=1.0)
    @test size(x) == (20, 20)
    @test sum(x) == 400
end

@testset "mask" begin
    A1 = [missing 1; 2 3]
    A2 = view([0 missing 1; 0 2 3], :, 2:3)
    ga1 = Raster(A1, (X, Y); missingval=missing)
    ga2 = Raster(A2, (X, Y); missingval=missing)
    @test all(mask(ga1; with=ga) .=== mask(ga2; with=ga) .=== [missing 1; 2 missing])
    ga2 = replace_missing(ga1 .* 1.0; missingval=NaN)
    @test all(mask(ga2; with=ga) .=== [NaN 1.0; 2.0 NaN])
    ga3 = replace_missing(ga1; missingval=-9999)
    mask!(ga3; with=ga)
    @test all(ga3 .=== [-9999 1; 2 -9999])
    dmask = mask(ga3; with=ga, filename="mask.tif")
    @test Rasters.isdisk(dmask)
    rm("mask.tif")
    stmask = mask(replace_missing(st, NaN); with=ga, filename="mask.tif")
    @test Rasters.isdisk(stmask)
    rm("mask_a.tif")
    rm("mask_b.tif")
    poly = polygon
    @testset "to polygon" begin
        for poly in (polygon, multi_polygon) 
            a1 = Raster(ones(X(-20:5; sampling=Intervals(Center())), Y(0:30; sampling=Intervals(Center()))))
            st1 = RasterStack(a1, a1)
            ser1 = RasterSeries([a1, a1], Ti(1:2))
            @test all(
                mask(a1; with=polygon) .===
                mask(st1; with=polygon)[:layer1] .===
                mask(ser1; with=polygon)[1]
            )
            # TODO: investigate this more for Points/Intervals
            # Exactly how do we define when boundary values are inside/outside a polygon
            @test sum(skipmissing(mask(a1; with=polygon, boundary=:inside))) == 19 * 19
            @test sum(skipmissing(mask(a1; with=polygon, boundary=:center))) == 20 * 20
            @test sum(skipmissing(mask(a1; with=polygon, boundary=:touches))) == 21 * 21
            mask(a1; with=polygon, boundary=:touches)
            mask(a1; with=polygon, boundary=:touches, shape=:line)
            mask(a1; with=polygon, boundary=:inside)
            mask(a1; with=polygon)
        end
    end
end

@testset "zonal" begin
    a = Raster((1:26) * (1:31)', (X(-20:5), Y(0:30)))
    zonal(sum, a; of=polygon) ==
        zonal(sum, a; of=[polygon, polygon])[1] ==
        zonal(sum, a; of=(geometry=polygon, x=:a, y=:b)) ==
        zonal(sum, a; of=[(geometry=polygon, x=:a, y=:b)])[1]
        zonal(sum, a; of=[(geometry=polygon, x=:a, y=:b)])[1] ==
        sum(skipmissing(mask(a; with=polygon)))
    @test zonal(sum, a; of=a) == 
        zonal(sum, a; of=dims(a)) ==
        zonal(sum, a; of=Extents.extent(a)) == 
        sum(a)

    b = a .* 2
    c = cat(a .* 3, a .* 3; dims=:newdim)
    st = RasterStack((; a, b, c))
    @test zonal(sum, st; of=polygon) == zonal(sum, st; of=[polygon])[1] ==
        zonal(sum, st; of=(geometry=polygon, x=:a, y=:b)) ==
        zonal(sum, st; of=[(geometry=polygon, x=:a, y=:b)])[1] ==
        zonal(sum, st; of=[(geometry=polygon, x=:a, y=:b)])[1] ==
        map(sum ∘ skipmissing, mask(st; with=polygon))  
    @test zonal(sum, st; of=st) == 
        zonal(sum, st; of=dims(st)) == 
        zonal(sum, st; of=Extents.extent(st)) == 
        sum(st)
end

@testset "classify" begin
    A1 = [missing 1; 2 3]
    ga1 = Raster(A1, (X, Y); missingval=missing)
    @test all(classify(ga1, 1=>99, 2=>88, 3=>77) .=== [missing 99; 88 77])
    @test all(classify(ga1, 1=>99, 2=>88, 3=>77; others=0) .=== [missing 99; 88 77])
    @test all(classify(ga1, 1=>99, 2=>88; others=0) .=== [missing 99; 88 0])
    A2 = [1.0 2.5; 3.0 4.0]
    ga2 = Raster(A2, (X, Y); missingval=missing)
    @test classify(ga2, (2, 3)=>:x, >(3)=>:y) == [1.0 :x; 3.0 :y]
    @test classify(ga2, (>=(1), <(2))=>:x, >=(3)=>:y) == [:x 2.5; :y :y]
    classify!(ga2, (1, 2.5)=>0.0, >=(3)=>-1.0; lower=(>), upper=(<=))
    @test ga2 == [1.0 0.0; -1.0 -1.0]
    classify!(ga2, [1 2.5 0.0; 2.5 4.0 -1.0]; lower=(>), upper=(<=))
    @test ga2 == [1.0 0.0; -1.0 -1.0]
    @test all(classify(ga1, [1 99; 2 88; 3 77]) .=== [missing 99; 88 77])
    @test_throws ArgumentError classify(ga1, [1, 2, 3])
end

@testset "points" begin
    ga = Raster(A, (X(9.0:1.0:10.0), Y(0.1:0.1:0.2)); missingval=missing)
    @test all(collect(points(ga; order=(Y, X))) .=== [missing (0.2, 9.0); (0.1, 10.0) missing])
    @test all(collect(points(ga; order=(X, Y))) .=== [missing (9.0, 0.2); (10.0, 0.1) missing])
    @test all(points(ga; order=(X, Y), ignore_missing=true) .===
              [(9.0, 0.1) (9.0, 0.2); (10.0, 0.1) (10.0, 0.2)])
end

# Idea for generic constructors
# Polygon(mod, values) = Polygon(Val{Symbol(mod)}(), values)
# Polygon(::Val{:ArchGDAL}, values) = ArchGDAL.createpolygon(values)

createpoint(args...) = ArchGDAL.createpoint(args...)
createfeature(x::Tuple{<:Any,<:Any}) = NamedTuple{(:geometry,:test)}(x)
createfeature(x::Tuple{<:Any,<:Any,<:Any}) = NamedTuple{(:geometry,:test,:test2)}(x)
createfeature(::Missing) = missing

@testset "extract" begin
    dimz = (X(9.0:1.0:10.0), Y(0.1:0.1:0.2))
    ga = Raster([1 2; 3 4], dimz; name=:test, missingval=missing)
    ga2 = Raster([5 6; 7 8], dimz; name=:test2, missingval=missing)
    st = RasterStack(ga, ga2)
    @testset "from Raster" begin
        # Tuple points
        @test all(extract(ga, [missing, (9.0, 0.1), (9.0, 0.2), (10.0, 0.3)]) .=== 
                  createfeature.([missing, ((9.0, 0.1), 1), ((9.0, 0.2), 2), ((10.0, 0.3), missing)]))
        # NamedTuple (reversed) points
        @test all(extract(ga, [missing, (Y=0.1, X=9.0), (Y=0.2, X=10.0), (Y=0.3, X=10.0)]) |> collect .=== 
                  createfeature.([missing, ((Y=0.1, X=9.0), 1), ((Y=0.2, X=10.0), 4), ((Y=0.3, X=10.0), missing)]))
        # Vector points
        @test all(extract(ga, [[9.0, 0.1], [10.0, 0.2]]) .== createfeature.([([9.0, 0.1], 1), ([10.0, 0.2], 4)]))
        # ArchGDAL equality is broken
        # @test all(extract(ga, ArchGDAL.createmultipoint([[0.1, 9.0], [0.2, 10.0], [0.3, 10.0]])) .==
                  # createfeature.([(createpoint(0.1, 9.0), 1), (createpoint(0.2, 10.0), 4), (createpoint(0.3, 10.0), missing)]))
        # Extract a polygon
        p = ArchGDAL.createpolygon([[[8.0, 0.0], [11.0, 0.0], [11.0, 0.4], [8.0, 0.0]]])
        @test all(extract(ga, p) .=== 
            createfeature.([((X=9.0, Y=0.1), 1), ((X=10.0, Y=0.1), 3), ((X=10.0, Y=0.2), 4)]))
    end
    @testset "from stack" begin
        @test all(extract(st, [missing, (9.0, 0.1), (10.0, 0.2), (10.0, 0.3)]) |> collect .===
                  createfeature.([missing, ((9.0, 0.1), 1, 5), ((10.0, 0.2), 4, 8), ((10.0, 0.3), missing, missing)]))
    end
end

@testset "trim, crop, extend" begin
    A = [missing missing missing
         missing 2.0     0.5
         missing 1.0     missing]

    r_fwd = Raster(A, (X(1.0:1.0:3.0), Y(1.0:1.0:3.0)); missingval=missing)
    r_revX = reverse(r_fwd; dims=X)
    r_revY = reverse(r_fwd; dims=Y)

    ga = r_revY
    for ga in (r_fwd, r_revX, r_revY)
        # Test with missing on all sides
        ga_r = rot180(ga)
        trimmed = trim(ga)
        trimmed_r = trim(ga_r)
        trimmed_d = trim(ga_r)
        if ga === r_fwd
            @test all(trimmed .=== [2.0 0.5; 1.0 missing])
            @test all(trimmed_r .=== [missing 1.0; 0.5 2.0])
        end
        cropped = crop(ga; to=trimmed)
        _, cropped2 = crop(trimmed, ga)
        cropped_r = crop(ga_r; to=trimmed_r)
        @test all(cropped .=== trimmed)
        @test all(cropped_r .=== trimmed_r)
        extended = extend(cropped, ga)[1]
        extended_r = extend(cropped_r; to=ga_r)
        extended_d = extend(cropped; to=ga, filename="extended.tif")
        @test all(extended .=== ga) 
        @test all(extended_r .=== ga_r)

        @testset "to polygons" begin
            A1 = Raster(zeros(X(-20:-5; sampling=Points()), Y(0:30; sampling=Points())))
            A2 = Raster(ones(X(-20:-5; sampling=Points()), Y(0:30; sampling=Points())))
            A1crop1 = crop(A1; to=polygon)
            A1crop2, A2crop = crop(A1, A2; to=polygon)
            size(A1crop1)
            @test size(A1crop1) == size(A1crop2) == size(A2crop) == (16, 21)
            @test bounds(A1crop1) == bounds(A1crop2) == bounds(A2crop) == ((-20, -5), (10, 30))
            A1extend1 = extend(A1; to=polygon)
            A1extend2, A2extend = extend(A1, A2; to=polygon)
            @test all(A1extend1 .=== A1extend2)
            @test size(A1extend1) == size(A1extend2) == size(A2extend) == (21, 31)
            @test bounds(A1extend1) == bounds(A1extend2) == bounds(A2extend) == ((-20.0, 0.0), (0, 30))
        end
        @testset "to featurecollection and table" begin
            A1 = Raster(zeros(X(-20:-5; sampling=Points()), Y(0:30; sampling=Points())))
            featurecollection = map(GeoInterface.getpoint(polygon), vals) do geometry, v
                (; geometry, val1=v, val2=2.0f0v)
            end
            fccrop = crop(A1; to=featurecollection)
            table = DataFrame(featurecollection)
            tablecrop = crop(A1; to=table)
            @test size(fccrop) == size(tablecrop) == (16, 21)
            @test bounds(fccrop) == bounds(tablecrop) == ((-20, -5), (10, 30))
        end
    end

end

@testset "mosaic" begin
    reg1 = Raster([0.1 0.2; 0.3 0.4], (X(2.0:1.0:3.0), Y(5.0:1.0:6.0)))
    reg2 = Raster([1.1 1.2; 1.3 1.4], (X(3.0:1.0:4.0), Y(6.0:1.0:7.0)))
    irreg1 = Raster([0.1 0.2; 0.3 0.4], (X([2.0, 3.0]), Y([5.0, 6.0])))
    irreg2 = Raster([1.1 1.2; 1.3 1.4], (X([3.0, 4.0]), Y([6.0, 7.0])))

    span_x1 = Explicit(vcat((1.5:1.0:2.5)', (2.5:1.0:3.5)'))
    span_x2 = Explicit(vcat((2.5:1.0:3.5)', (3.5:1.0:4.5)'))
    exp1 = Raster([0.1 0.2; 0.3 0.4], (X(Sampled([2.0, 3.0]; span=span_x1)), Y([5.0, 6.0])))
    exp2 = Raster([1.1 1.2; 1.3 1.4], (X(Sampled([3.0, 4.0]; span=span_x2)), Y([6.0, 7.0])))
    @test val(span(mosaic(first, exp1, exp2), X)) == [1.5 2.5 3.5; 2.5 3.5 4.5]
    @test all(mosaic(first, [reg1, reg2]) .=== 
              mosaic(first, irreg1, irreg2) .===
              mosaic(first, (irreg1, irreg2)) .=== 
              [0.1 0.2 missing; 
               0.3 0.4 1.2; 
               missing 1.3 1.4])
    @test all(mosaic(last, reg1, reg2) .===
              mosaic(last, irreg1, irreg2) .===
              mosaic(last, exp1, exp2) .=== [0.1 0.2 missing; 
                                             0.3 1.1 1.2; 
                                             missing 1.3 1.4])

    @test all(mosaic(first, [reverse(reg2; dims=Y), reverse(reg1; dims=Y)]) .=== 
              [missing 0.2 0.1; 
               1.2 1.1 0.3; 
               1.4 1.3 missing]
             )

    # 3 dimensions
    A1 = Raster(ones(2, 2, 2), (X(2.0:-1.0:1.0), Y(5.0:1.0:6.0), Ti(DateTime(2001):Year(1):DateTime(2002))))
    A2 = Raster(zeros(2, 2, 2), (X(3.0:-1.0:2.0), Y(4.0:1.0:5.0), Ti(DateTime(2002):Year(1):DateTime(2003))))
    @test all(mosaic(mean, A1, A2) |> parent .=== 
              first(mosaic(mean, RasterStack(A1), RasterStack(A2))) .===
              cat([missing missing missing
                   missing 1.0     1.0
                   missing 1.0     1.0    ],
                   [0.0     0.0 missing
                   0.0     0.5     1.0   
                   missing 1.0     1.0    ],
                   [0.0     0.0     missing
                   0.0     0.0     missing    
                   missing missing missing], dims=3))
end

@testset "rasterize and rasterize!" begin
    A1 = Raster(zeros(X(-20:5; sampling=Intervals()), Y(0:30; sampling=Intervals())))
    A2 = Raster(zeros(Y(0:30; sampling=Intervals()), X(-20:5; sampling=Intervals())))
    st = RasterStack((A1, copy(A1)))

    @testset "all geoms work as :point" begin
        A = A1
        geom = pointvec
        for A in (A1, A2), geom in (pointvec, pointfc, multi_point, linestring, multi_linestring, linearring, polygon, multi_polygon)
            A .= 0
            rasterize!(A, geom; shape=:point, fill=1)
            st.layer1 .= st.layer2 .= 0
            rasterize!(st, geom; shape=:point, fill=(layer1=1, layer2=2))
            @test sum(st[:layer1]) == 4
            @test sum(st[:layer2]) == 8
            st[:layer1] .= st[:layer2] .= 0
            @test_nowarn rasterize!(A, geom; shape=:point, fill=1)
        end
    end

    @testset "all line and polygon geoms work as :line" begin
        A = A1
        geom = linestring
        for A in (A1, A2), geom in (linestring, multi_linestring, linearring, polygon, multi_polygon)
            rasterize!(A, geom; shape=:line, fill=1)
            @test sum(A) == 20 + 20 + 20 + 20
            A .= 0
        end
        @testset ":line is detected for line geometries" begin
            for A in (A1, A2), geom in (linestring, multi_linestring)
                rasterize!(A, geom; fill=1)
                @test sum(A) == 20 + 20 + 20 + 20
                A .= 0
            end
        end
    end

    @testset "polygon geoms work as :polygon" begin
        A = A1
        poly = polygon
        for A in (A1, A2), poly in (polygon, multi_polygon)
            ra = rasterize(poly; to=A, missingval=0, shape=:polygon, fill=1, boundary=:center)
            ra_res = rasterize(poly; res=map(step, span(A)), missingval=0, shape=:polygon, fill=1, boundary=:center)
            @test sum(ra) == sum(ra_res) === 20 * 20
            ra = rasterize(poly; to=A, shape=:polygon, fill=1, boundary=:touches)
            @test sum(skipmissing(ra)) === 21 * 21
            rasterize!(A, poly; shape=:polygon, fill=1, boundary=:inside)
            @test sum(A) === 19.0 * 19.0
            A .= 0

            @testset "polygon is detected for polygon geometries" begin
                A = Raster(zeros(X(-20:5; sampling=Intervals()), Y(0:30; sampling=Intervals())))
                R = rasterize(poly; to=A, fill=1)
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
                rst = rasterize(feature; to=A, fill=(:val1, :val2))
                @test keys(rst) == (:val1, :val2)
                @test dims(rst) == dims(A)
                @test map(sum ∘ skipmissing, rst) === (val1=4, val2=8.0f0)
            end
            @testset "Symbol fill makes an array" begin
                ra = rasterize(feature; to=A, fill=:val1)
                @test ra isa Raster{Union{Missing,Int64}}
                @test name(ra) == :val1
            end
        end

        @testset "feature collection, table from fill of Symbol keys" begin
            data = pointfc
            pointfc
            for data in (pointfc, DataFrame(pointfc))
                @testset "NTuple of Symbol fill makes an stack" begin
                    rst = rasterize(sum, data; to=A, fill=(:val1, :val2))
                    @test keys(rst) == (:val1, :val2)
                    @test map(sum ∘ skipmissing, rst) === (val1=14, val2=28.0f0)
                    @test_throws ArgumentError rasterize(data; to=A, fill=(:val1, :not_a_column))
                end
                @testset "Symbol fill makes an array" begin
                    ra = rasterize(data; to=A, fill=:val1)
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
        pointvec2 = map(p -> (p[1] + 10, p[2] + 10), pointvec)
        pointvec3 = map(p -> (p[1] + 20, p[2] + 20), pointvec)
        pointvec4 = map(p -> (p[1] + 30, p[2] + 30), pointvec)
        polygon = ArchGDAL.createpolygon(pointvec)
        polygons = ArchGDAL.createpolygon.([[pointvec], [pointvec2], [pointvec3], [pointvec4]])
        # With fill of these are all the same thing
        reduced_raster_last_center = rasterize(last, polygons; res=5, fill=1, boundary=:center)
        reduced_raster_first_center = rasterize(first, polygons; res=5, fill=1, boundary=:center)
        reduced_raster_mean_center = rasterize(mean, polygons; res=5, fill=1, boundary=:center)
        reduced_raster_median_center = rasterize(median, polygons; res=5, fill=1, boundary=:center)
        # using Plots
        # plot(reduced_raster_mean_center)
        # plot(reduced_raster_last_center)
        # plot(reduced_raster_median_center)
        @test sum(skipmissing(reduced_raster_last_center)) ==
              sum(skipmissing(reduced_raster_first_center)) ==
              sum(skipmissing(reduced_raster_median_center)) ==
              sum(skipmissing(reduced_raster_mean_center)) == 16 * 4 - 12
        # The outlines of these plots should exactly mactch,  
        # Plots.plot(raster; clims=(0, 3))
        # Plots.plot!(polygons; opacity=0.3, fillcolor=:black)
        
        
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
end

@testset "coverage" begin
    @time covsum = Rasters.coverage(shphandle.shapes; mode=sum, res=1, scale=10);
    @time covunion = Rasters.coverage(shphandle.shapes; mode=union, res=1, scale=10);
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

@testset "resample" begin
    raster_path = maybedownload("https://download.osgeo.org/geotiff/samples/gdal_eg/cea.tif")

    output_res = 0.0027
    output_crs = EPSG(4326)
    resample_method = "near"

    ## Resample cea.tif manually with ArchGDAL
    wkt = convert(String, convert(WellKnownText, output_crs))
    AG_output = ArchGDAL.read(raster_path) do dataset
        ArchGDAL.gdalwarp([dataset], ["-t_srs", "$(wkt)",
                                "-tr", "$(output_res)", "$(output_res)",
                                "-r", "$(resample_method)"]) do warped
            ArchGDAL.read(ArchGDAL.getband(warped, 1))
        end
    end

    # Resample cea.tif using resample
    cea = Raster(raster_path)
    raster_output = resample(cea, output_res; crs=output_crs, method=resample_method)
    disk_output = resample(cea, output_res; crs=output_crs, method=resample_method, filename="resample.tif")

    cea_permuted = permutedims(Raster(raster_path), (Y, X, Band))
    permuted_output = resample(cea_permuted, output_res; crs=output_crs, method=resample_method)

    # Compare ArchGDAL, resample and permuted resample 
    @test AG_output ==
        raster_output[Band(1)] ==
        disk_output[Band(1)] ==
        permutedims(permuted_output, (X, Y, Band))[Band(1)]
    @test abs(step(dims(raster_output, Y))) ≈
        abs(step(dims(raster_output, X))) ≈ 
        abs(step(dims(disk_output, X))) ≈ 
        abs(step(dims(permuted_output, X))) ≈ output_res

    rm("resample.tif")

    @testset "snapped size and dim index match" begin
        snaptarget = raster_output
        snapped = resample(cea; to=snaptarget)
        disk_snapped = resample(cea; to=snaptarget, filename="raster.tif")
        @test size(snapped) == size(disk_snapped) == size(snaptarget)
        @test isapprox(index(snaptarget, Y), index(snapped, Y))
        @test isapprox(index(snaptarget, X), index(snapped, X))
        @test isapprox(index(snaptarget, Y), index(disk_snapped, Y))
        @test isapprox(index(snaptarget, X), index(disk_snapped, X))
    end
end
