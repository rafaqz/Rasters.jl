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

test_shape_dir = realpath(joinpath(dirname(pathof(Shapefile)), "..", "test", "shapelib_testcases"))
shp_paths = filter(x -> occursin("shp", x), readdir(test_shape_dir; join=true))
shppath = shp_paths[1]
shphandle = Shapefile.Handle(shppath)

ga99 = replace_missing(ga, -9999)
gaNaN = replace_missing(ga, NaN32)
gaMi = replace_missing(ga)

@testset "replace_missing" begin
    @test eltype(ga99) == Float32
    @test eltype(gaNaN) == Float32
    @test eltype(gaMi) == Union{Float32,Missing}
    @test eltype(replace_missing(ga, 0.0)) == Float64
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
    @test parent(boolmask(ga)) isa BitMatrix
    @test boolmask(ga99) == [false true; true false]
    @test boolmask(gaNaN) == [false true; true false]
    @test dims(boolmask(ga)) == (X(NoLookup(Base.OneTo(2))), Y(NoLookup(Base.OneTo(2))))
    x = boolmask(polygon; res=1.0) 
    @test x == trues(X(Projected(-20:1.0:-1.0; crs=nothing)), Y(Projected(10.0:1.0:29.0; crs=nothing)))
    @test parent(x) isa BitMatrix
    # With a :geometry axis
    x = boolmask([polygon, polygon]; collapse=false, res=1.0)
    @test eltype(x) == Bool
    @test size(x) == (20, 20, 2)
    @test sum(x) == 800
    @test parent(x) isa BitArray{3}
    x = boolmask([polygon, polygon]; collapse=true, res=1.0)
    @test size(x) == (20, 20)
    @test sum(x) == 400
    @test parent(x) isa BitMatrix
end

@testset "missingmask" begin
    @test all(missingmask(ga) .=== [missing true; true missing])
    @test all(missingmask(ga99) .=== [missing true; true missing])
    @test all(missingmask(gaNaN) .=== [missing true; true missing])
    @test dims(missingmask(ga)) == (X(NoLookup(Base.OneTo(2))), Y(NoLookup(Base.OneTo(2))))
    @test missingmask(polygon; res=1.0) == fill!(Raster{Union{Missing,Bool}}(undef, X(Projected(-20:1.0:-1.0; crs=nothing)), Y(Projected(10.0:1.0:29.0; crs=nothing))), true)
    x = missingmask([polygon, polygon]; collapse=false, res=1.0)
    @test eltype(x) == Union{Bool,Missing}
    @test size(x) == (20, 20, 2)
    @test sum(x) == 800
    @test parent(x) isa Array{Union{Missing,Bool},3}
    x = missingmask([polygon, polygon]; collapse=true, res=1.0)
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

@testset "mask_replace_missing" begin
    # Floating point rasters
    a = Raster([1.0 0.0; 1.0 1.0], dims=(X, Y), missingval=0)
    b = Raster([1.0 1.0; 1.0 0.0], dims=(X, Y), missingval=0)

    # Integer rasters
    c = Raster([1 0; 1 1], dims=(X, Y), missingval=0)
    d = Raster([1 1; 1 0], dims=(X, Y), missingval=0)

    # Test that missingval is replaced in source mask (Floats)
    @test isequal(mask(a, with=b, missingval=3.14), [1.0 3.14; 1.0 3.14]) # Test missingval = 3.14
    @test isequal(mask(a, with=b, missingval=missing), [1.0 missing; 1.0 missing]) # Test missingval = missing
    @test isequal(mask(a, with=b, missingval=NaN), [1.0 NaN; 1.0 NaN]) # Test missingval = NaN
    @test isequal(mask(a, with=b, missingval=NaN32), [1.0 NaN; 1.0 NaN]) # Test convert NaN32 to NaN
    @test isequal(mask(a, with=b, missingval=Inf), [1.0 Inf; 1.0 Inf]) # Test missingval = Inf
    @test_throws MethodError mask(a, with=b, missingval=nothing)

    # Test that missingval is replaced in source mask (Ints)
    @test isequal(mask(c, with=d, missingval=missing), [1 missing; 1 missing]) # Test missingval = missing
    @test isequal(mask(c, with=d, missingval=-1.0), [1 -1; 1 -1])
    @test_throws MethodError mask(c, with=d, missingval=nothing)
    @test_throws InexactError mask(c, with=d, missingval=NaN)
    @test_throws InexactError mask(c, with=d, missingval=3.14)
    @test_throws InexactError mask(c, with=d, missingval=Inf)

    # Test Type Stability
    @test eltype(mask(a, with=b, missingval=0)) == Float64
    @test eltype(mask(a, with=b, missingval=-1)) == Float64
    @test eltype(mask(a, with=b, missingval=Inf32)) == Float64
    @test eltype(mask(Float32.(a), with=b, missingval=Inf)) == Float32
    @test eltype(mask(Float32.(a), with=b, missingval=NaN)) == Float32
    @test eltype(mask(Float32.(a), with=b, missingval=0.0)) == Float32
    @test eltype(mask(Float32.(a), with=b, missingval=0)) == Float32
    @test eltype(mask(Float32.(a), with=b, missingval=-1)) == Float32
    @test eltype(mask(c, with=d, missingval=-1.0)) == Int64
    @test eltype(mask(c, with=d, missingval=0.0f0)) == Int64
    @test eltype(mask(c, with=Float64.(d), missingval=-1.0)) == Int64
    @test eltype(mask(c, with=Float64.(d), missingval=0.0f0)) == Int64

    # Test mask!
    @test_throws MethodError mask!(a, with=b, missingval=missing)
    @test isequal(mask!(deepcopy(a), with=b, missingval=3.14), [1.0 3.14; 1.0 3.14]) # Test missingval = 3.14
    @test isequal(mask!(deepcopy(a), with=b, missingval=NaN), [1.0 NaN; 1.0 NaN]) # Test missingval = NaN
    @test isequal(mask!(deepcopy(a), with=b, missingval=NaN32), [1.0 NaN; 1.0 NaN]) # Test convert NaN32 to NaN
    @test isequal(mask!(deepcopy(a), with=b, missingval=Inf), [1.0 Inf; 1.0 Inf]) # Test missingval = Inf
    @test isequal(mask(deepcopy(c), with=d, missingval=-1.0), [1 -1; 1 -1])
    @test_throws MethodError mask!(deepcopy(a), with=b, missingval=missing)
    @test_throws ArgumentError mask!(deepcopy(c), with=d, missingval=nothing)
    @test_throws InexactError mask!(deepcopy(c), with=d, missingval=NaN)
    @test_throws InexactError mask!(deepcopy(c), with=d, missingval=3.14)
    @test_throws InexactError mask!(deepcopy(c), with=d, missingval=Inf)
end

@testset "zonal" begin
    a = Raster((1:26) * (1:31)', (X(-20:5), Y(0:30)))
    pointvec_empty = [(-100.0, 0.0), (-100.0, 0.0), (-100.0, 0.0), (-100.0, 0.0), (-100.0, 0.0)]
    polygon_empty = ArchGDAL.createpolygon(pointvec_empty)
    zonal(sum, a; of=polygon) ==
        zonal(sum, a; of=[polygon, polygon])[1] ==
        zonal(sum, a; of=[polygon, polygon_empty])[1] ==
        zonal(sum, a; of=(geometry=polygon, x=:a, y=:b)) ==
        zonal(sum, a; of=[(geometry=polygon, x=:a, y=:b)])[1] ==
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

@testset "resample" begin
    raster_path = maybedownload("https://download.osgeo.org/geotiff/samples/gdal_eg/cea.tif")

    output_res = 0.0027
    output_crs = EPSG(4326)
    method = "near"

    ## Resample cea.tif manually with ArchGDAL
    wkt = convert(String, convert(WellKnownText, output_crs))
    AG_output = ArchGDAL.read(raster_path) do dataset
        ArchGDAL.gdalwarp([dataset], ["-t_srs", "$(wkt)",
                                "-tr", "$(output_res)", "$(output_res)",
                                "-r", "$(method)"]) do warped
            ArchGDAL.read(ArchGDAL.getband(warped, 1))
        end
    end

    # Resample cea.tif using resample
    cea = Raster(raster_path; missingval=0x00)
    raster_output = 
    resample(cea; res=output_res, crs=output_crs, method)
    disk_output = resample(cea; res=output_res, crs=output_crs, method, filename="resample.tif")

    cea_permuted = permutedims(Raster(raster_path), (Y, X))
    permuted_output = resample(cea_permuted, output_res; crs=output_crs, method)

    # Compare ArchGDAL, resample and permuted resample 
    @test AG_output ==
        raster_output[Band(1)] ==
        disk_output[Band(1)] ==
        permutedims(permuted_output, (X, Y))
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

    @testset "`method` only does nothing" begin
        resampled = resample(cea; method)
        @test crs(cea) == crs(resampled)
        @test cea == resampled
        # There is some floating point error here after Rasters -> GDAL -> Rasterss...
        # Should we correct it by detecting almost identical extent and using the original?
        @test_broken extent(cea) = extent(resampled)
    end

    @testset "only `res` kw changes the array size predictably" begin
        res = step(span(cea, X)) / 2
        resampled = resample(cea; res)
        @test crs(cea) == crs(resampled)
        @test size(dims(resampled, (X, Y))) == size(dims(cea, (X, Y))) .* 2
        # GDAL fp error see above
        @test_broken extent(cea) = extent(resampled)
        resampled = resample(cea; res=(res, 2res))
        @test size(dims(resampled, (X, Y))) == (size(cea, X) .* 2, size(cea, Y))
        resampled = resample(cea; res=(X(2res), Y(res)))
        @test size(dims(resampled, (X, Y))) == (size(cea, X), size(cea, Y) * 2)
    end

    @testset "only `size` kw sets the size" begin
        res = step(span(cea, X)) / 2
        resampled = resample(cea; size=(100, 200))
        @test crs(cea) == crs(resampled)
        @test size(dims(resampled, (X, Y))) == size(resampled[:, :, 1]) == (100, 200)
        resampled = resample(cea; size=(X(99), Y(111)))
        @test crs(cea) == crs(resampled)
        @test size(dims(resampled, (X, Y))) == size(resampled[:, :, 1]) == (99, 111)
        resampled = resample(cea; size=100)
        @test size(dims(resampled, (X, Y))) == size(resampled[:, :, 1]) == (100, 100)
    end

    @testset "Extent `to` can resize arbitrarily" begin
        to = Extent(X=(-10000.0, 2000.0), Y=(4.2e6, 4.3e6))
        resampled = resample(cea; to)
        @test all(map((bs...,) -> all(map(≈, bs...)), extent(resampled), to))
        @test crs(cea) == crs(resampled)
        # Most of the area is extended, and filled with missingval
        @test sum(x -> x === missingval(resampled), resampled) > length(resampled) / 2
    end

    @testset "Geometries work as `to`" begin
        geom = GeoInterface.Polygon([[(0.0, 4e6), (-1e5, 4.4e6), (-1e5, 4.2e6)]])
        resampled = resample(cea; to=geom)
        @test map(extent(resampled[Band(1)]), GeoInterface.extent(geom)) do bs1, bs2 
            map(≈, bs2, bs2) |> all 
        end |> all
        @test crs(cea) == crs(resampled)
        # Most of the area is extended, and filled with missingval
        @test sum(x -> x === missingval(resampled), resampled) > length(resampled) / 2
    end

    @testset "only `crs` kw changes the array size" begin
        resampled = resample(cea; crs=EPSG(3857), method)
        @test size(dims(resampled, (X, Y))) !== size(dims(cea, (X, Y)))
        @test crs(resampled) == EPSG(3857)
    end

    @testset "no existing crs warns" begin
        nocrs = setcrs(cea, nothing)
        @test crs(nocrs) == nothing
        @test_warn "does not have crs" resample(nocrs; crs=output_crs, method)
    end

    @testset "resample eltype propagates" begin
        r = Raster(rand(UInt8, X(1:10), Y(1:10)))
        r1 = resample(r; to=r)
        @test eltype(r1) == UInt8
    end
end
