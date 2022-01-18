using Rasters, Test, ArchGDAL, Dates, Statistics, GeoInterface
using Rasters.LookupArrays, Rasters.Dimensions 
using GeoInterface: Point, Polygon, LineString

include(joinpath(dirname(pathof(Rasters)), "../test/test_utils.jl"))

A = [missing 7.0f0; 2.0f0 missing]
B = [1.0 0.4; 2.0 missing]
ga = Raster(A, (X, Y); missingval=missing)
ga99 = replace_missing(ga, -9999)
gaNaN = replace_missing(ga, NaN32)
gaMi = replace_missing(ga)
st = RasterStack((a=A, b=B), (X, Y); missingval=(a=missing,b=missing))

polygon = [[-20.0, 30.0],
           [-20.0, 10.0],
           [0.0, 10.0],
           [0.0, 30.0],
           [-20.0, 30.0]]
vals = [1, 2, 3, 4, 5]

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
    @test Rasters.isdisk(dNaN)
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
end

@testset "missingmask" begin
    @test all(missingmask(ga) .=== [missing true; true missing])
    @test all(missingmask(ga99) .=== [missing true; true missing])
    @test all(missingmask(gaNaN) .=== [missing true; true missing])
    @test dims(missingmask(ga)) == (X(NoLookup(Base.OneTo(2))), Y(NoLookup(Base.OneTo(2))))
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
    @testset "to polygon" begin
        a1 = Raster(ones(X(-20:5), Y(0:30)))
        st1 = RasterStack(a1, a1)
        ser1 = RasterSeries([a1, a1], Ti(1:2))
        @test all(mask(a1; with=polygon) .===
            mask(st1; with=polygon)[:layer1] .===
            mask(ser1; with=polygon)[1]
        )
        # TODO: investigate this more for Points/Intervals
        # Exactly how do we define when boundary values are inside/outside a polygon
        @test sum(skipmissing(mask(a1; with=polygon, boundary=:inside))) == 19 * 19
        @test sum(skipmissing(mask(a1; with=polygon, boundary=:center))) == 20 * 20
        @test sum(skipmissing(mask(a1; with=polygon, boundary=:touches))) == 21 * 21
    end
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
    @test all(
              collect(
                      points(ga; order=(Y, X)))
              .=== [missing (0.2, 9.0); (0.1, 10.0) missing])
    @test all(collect(points(ga; order=(X, Y))) .=== [missing (9.0, 0.2); (10.0, 0.1) missing])
    @test all(points(ga; order=(X, Y), ignore_missing=true) .===
              [(9.0, 0.1) (9.0, 0.2); (10.0, 0.1) (10.0, 0.2)])
end

@testset "extract" begin
    A1 = [1 2; 3 4]
    A2 = [5 6; 7 8]
    dimz = (X(9.0:1.0:10.0), Y(0.1:0.1:0.2))
    ga = Raster(A1, dimz; name=:test, missingval=missing)
    ga2 = Raster(A2, dimz; name=:test2, missingval=missing)
    st = RasterStack(ga, ga2)
    @testset "from Raster" begin
        @test all(extract(ga, [missing, [9.0, 0.1], [10.0, 0.2], [10.0, 0.3]]) .=== 
                  [missing, (X=9.0, Y=0.1, test=1), (X=10.0, Y=0.2, test=4), (X=10.0, Y=0.3, test=missing)])
        @test all(extract(ga, ([9.0, 10.0, 10.0], [0.1, 0.2, 0.3])) .=== 
                  [(X=9.0, Y=0.1, test=1), (X=10.0, Y=0.2, test=4), (X=10.0, Y=0.3, test=missing)])
        @test all(extract(ga, Polygon([[9.0, 0.1], [10.0, 0.2], [10.0, 0.3]])) .=== 
                  [(X=9.0, Y=0.1, test=1), (X=10.0, Y=0.2, test=4), (X=10.0, Y=0.3, test=missing)])
        @test all(extract(ga, [(0.1, 9.0), (0.2, 10.0), (0.3, 10.0), missing]; order=(Y, X)) .=== 
                  [(Y=0.1, X=9.0, test=1), (Y=0.2, X=10.0, test=4), (Y=0.3, X=10.0, test=missing), missing])
        @test all(extract(ga, ([0.1, 0.2, 0.3], [9.0, 10.0, 10.0]); order=(Y, X)) .=== 
                  [(Y=0.1, X=9.0, test=1), (Y=0.2, X=10.0, test=4), (Y=0.3, X=10.0, test=missing)])
        @test all(extract(ga, Polygon([[0.1, 9.0], [0.2, 10.0], [0.3, 10.0]]); order=(Y, X)) .=== 
                  [(Y=0.1, X=9.0, test=1), (Y=0.2, X=10.0, test=4), (Y=0.3, X=10.0, test=missing)])
        @test all(extract(ga, [(missing, 9.0), (0.2, 10.0), (0.3, 10.0)]; order=(Y, X)) .=== 
                  [(Y=missing, X=missing, test=missing), (Y=0.2, X=10.0, test=4), (Y=0.3, X=10.0, test=missing)])
        @test all(extract(ga, [[missing, 9.0], [0.2, 10.0], [0.3, 10.0]]; order=(Y, X)) .=== 
                  [(Y=missing, X=missing, test=missing), (Y=0.2, X=10.0, test=4), (Y=0.3, X=10.0, test=missing)])
        @test_throws ArgumentError extract(ga, [missing, [9.0, 0.1], [10.0, 0.2], [10.0, 0.3]]; order=())
        @test_throws ArgumentError extract(ga, [missing, [9.0, 0.1], [10.0, 0.2], [10.0, 0.3]]; order=(Z, X))
    end
    @testset "from stack" begin
        @test all(extract(st, [missing, [9.0, 0.1], [10.0, 0.2], [10.0, 0.3]]) .=== 
              [missing, (X=9.0, Y=0.1, test=1, test2=5), (X=10.0, Y=0.2, test=4, test2=8), (X=10.0, Y=0.3, test=missing, test2=missing)])
    end
    @testset "Tables.jl compatible" begin
        @test all(extract(ga, [(X=9.0, Y=0.1), (X=10.0, Y=0.2), (X=10.0, Y=0.3)]) .=== 
              [(X=9.0, Y=0.1, test=1), (X=10.0, Y=0.2, test=4), (X=10.0, Y=0.3, test=missing)])
        @test all(extract(ga, [(X=9.0, Y=0.1), (X=10.0, Y=0.2), (X=10.0, Y=0.3)]; order=(Y=>:Y, X=>:X)) .=== 
              [(Y=0.1, X=9.0, test=1), (Y=0.2, X=10.0, test=4), (Y=0.3, X=10.0, test=missing)])
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
        @testset "to table" begin
            A1 = Raster(zeros(X(-20:-5; sampling=Points()), Y(0:30; sampling=Points())))
            table = map(polygon, vals) do p, v
                (x=p[1], y=p[2], val1=v, val2=2.0f0v)
            end
            A1crop1 = crop(A1; to=table, order=(X=>:x, Y=>:y))
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

@testset "inpolygon" begin
    gi_polygon = Polygon([polygon])
    poly = polygon
    for poly in (polygon, gi_polygon)
        @test inpolygon((-10.0, 20.0), poly) == true
        @test inpolygon((-19.0, 29.0), poly) == true
        @test inpolygon((-30.0, 20.0), poly) == false
        @test inpolygon([(-10.0, 20.0), (-30.0, 40.0)], poly) == [true, false]
        @test inpolygon(Point([-20.0, 50.0]), poly) == false
        @test inpolygon(LineString([[-10.0, 20.0], [-30.0, 40.0]]), poly) == [true, false]
        @test inpolygon(Polygon([[[-10.0, 20.0], [-30.0, 40.0]]]), poly) == [true, false]
    end
end

@testset "rasterize" begin
    gi_polygon = Polygon([polygon])
    rev_polygon = reverse.(polygon)
    A1 = Raster(zeros(X(-20:5; sampling=Intervals()), Y(0:30; sampling=Intervals())))
    A2 = Raster(zeros(Y(0:30; sampling=Intervals()), X(-20:5; sampling=Intervals())))
    st = RasterStack((A1, A1))

    A = A1
    poly = polygon
    poly = gi_polygon
    ord = (X, Y)
    for A in (A1, A2), (ord, poly) in (((X, Y), polygon), ((X, Y), gi_polygon), ((Y, X), rev_polygon))
        A .= 0
        rasterize!(A, poly, vals; order=ord)
        @test sum(A) == 14 # The last value overwrites the first
        A .= 0
        rasterize!(A, poly; shape=:point, fill=1, order=ord)
        @test sum(A) == 4
        A .= 0
        rasterize!(A, poly; shape=:line, fill=1, order=ord)
        @test sum(A) == 20 + 20 + 20 + 20
        A .= 0
        rasterize!(A, poly; shape=:polygon, fill=1, boundary=:center, order=ord)
        @test sum(A) == 20 * 20
        A .= 0
        rasterize!(A, poly; shape=:polygon, fill=1, boundary=:touches, order=ord)
        @test sum(A) == 21 * 21
        A .= 0
        rasterize!(A, poly; shape=:polygon, fill=1, boundary=:inside, order=ord)
        @test sum(A) == 19 * 19
        A = Raster(zeros(X(-20:5; sampling=Intervals()), Y(0:30; sampling=Intervals())))
        R = rasterize(poly, vals; to=A, order=ord)
        @test sum(skipmissing(R)) == 14 # The last value overwrites the first
        R = rasterize(poly; fill=1, to=A, order=ord)
        @test sum(skipmissing(R)) == 20 * 20
        @test_throws ArgumentError rasterize(poly, vals; to=A, order=())
        @test_throws ArgumentError rasterize!(A, poly; shape=:notashape, fill=1, order=ord)
        @test_throws ArgumentError rasterize!(A, poly; shape=:polygon, fill=1, boundary=:notaboundary, order=ord)

        st[:layer1] .= st[:layer2] .= 0
        rasterize!(st, poly; shape=:point, fill=1, order=ord)
        @test sum(st[:layer1]) == sum(st[:layer2]) == 4
        st[:layer1] .= st[:layer2] .= 0
        ra1 = rasterize(poly, vals; to=st, order=ord)
        @test name(ra1) == :layer1
        @test sum(skipmissing(ra1)) == 14 # The last value overwrites the first
        stvals = [(1, 2), (3, 4), (5, 6), (7, 8), (9, 10)]
        st2 = rasterize(poly, stvals; to=st, order=ord)
        @test sum(skipmissing(st2[:layer1])) == 24 # The last value overwrites the first
        @test sum(skipmissing(st2[:layer2])) == 28
    end

    @testset "table" begin
        table = map(polygon, vals) do p, v
            (x=p[1], y=p[2], val1=v, val2=2.0f0v)
        end
        @test sum(skipmissing(rasterize(table; to=A, order=(X=>:x, Y=>:y), name=:val1))) == 14
        rst = rasterize(table; to=A, order=(X=>:x, Y=>:y), name=(:val1, :val2))
        @test map(sum ∘ skipmissing, rst) === (val1=14, val2=28.0f0)
        R = rasterize(table; to=A, order=(X=>:x, Y=>:y)) 
        @test sum(skipmissing(R[:val1])) === 14
        @test sum(skipmissing(R[:val2])) === 28.0f0
        @test keys(R) == (:val1, :val2)
        @test dims(R) == dims(A)
        @test_throws ArgumentError rasterize(table; to=A, order=(), name=:val1)
        @test_throws ArgumentError rasterize(table; to=A, order=(Z,), name=:val1)
        @test_throws ArgumentError rasterize(table; to=A, name=:val1)
    end

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

    ## Resample cea.tif using resample
    cea = read(geoarray(raster_path))
    raster_output = resample(cea, output_res; crs=output_crs, method=resample_method)
    disk_output = resample(cea, output_res; crs=output_crs, method=resample_method, filename="resample.tif")

    cea_permuted = permutedims(read(geoarray(raster_path)), (Y, X, Band))
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
