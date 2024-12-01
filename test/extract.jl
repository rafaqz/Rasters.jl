using Rasters, Test, DataFrames, Extents
import GeoInterface as GI
using Rasters.Lookups, Rasters.Dimensions 

include(joinpath(dirname(pathof(Rasters)), "../test/test_utils.jl"))

dimz = (X(9.0:1.0:10.0), Y(0.1:0.1:0.2))
rast = Raster(Union{Int,Missing}[1 2; 3 4], dimz; name=:test, missingval=missing)
rast2 = Raster([5 6; 7 8], dimz; name=:test2, missingval=5)
rast_m = Raster([1 2; 3 missing], dimz; name=:test, missingval=missing)
st = RasterStack(rast, rast2)

points = [missing, (9.0, 0.1), (9.0, 0.2), (10.0, 0.3), (10.0, 0.2)]
poly = GI.Polygon([[(8.0, 0.0), (11.0, 0.0), (11.0, 0.4), (8.0, 0.0)]])
linestring = GI.LineString([(8.0, 0.0), (11.0, 0.0), (11.0, 0.4)])
line = GI.Line([(8.0, 0.0), (12.0, 4.0)])
table = (geometry=points, foo=zeros(4))

@testset "Points" begin
    @testset "From Raster" begin
        # Tuple points
        ex = extract(rast, points)
        T = @NamedTuple{geometry::Union{Missing,Tuple{Float64,Float64}},test::Union{Missing,Int64}}
        @test eltype(ex) == T
        @test all(ex .=== T[
            (geometry = missing, test = missing)
            (geometry = (9.0, 0.1), test=1)
            (geometry = (9.0, 0.2), test=2)
            (geometry = (10.0, 0.3), test=missing)
            (geometry = (10.0, 0.2), test=4)
        ])
        ex = extract(rast_m, points; skipmissing=true)
        T = @NamedTuple{geometry::Tuple{Float64, Float64}, test::Int64}
        @test eltype(ex) == T
        @test all(ex .=== T[(geometry = (9.0, 0.1), test = 1), (geometry = (9.0, 0.2), test = 2)])
        ex = extract(rast_m, points; skipmissing=true, geometry=false)
        T = @NamedTuple{test::Int64}
        @test eltype(ex) == T
        @test all(ex .=== T[(test = 1,), (test = 2,)])
        @test all(extract(rast_m, points; skipmissing=true, geometry=false, index=true) .=== [
            (index = (1, 1), test = 1,)
            (index = (1, 2), test = 2,)
        ])
        # NamedTuple (reversed) points - tests a Table that iterates over points
        T = @NamedTuple{geometry::Union{@NamedTuple{Y::Float64,X::Float64}},test::Union{Missing,Int64}}
        @test all(extract(rast, [(Y=0.1, X=9.0), (Y=0.2, X=10.0), (Y=0.3, X=10.0)]) .=== T[
            (geometry = (Y = 0.1, X = 9.0), test = 1)
            (geometry = (Y = 0.2, X = 10.0), test = 4)
            (geometry = (Y = 0.3, X = 10.0), test = missing)
        ])
        # Vector points
        @test all(extract(rast, [[9.0, 0.1], [10.0, 0.2]]) .== [
            (geometry = [9.0, 0.1], test = 1)
            (geometry = [10.0, 0.2], test = 4)
        ])
    end

    @testset "From RasterStack" begin
        T = @NamedTuple{geometry::Union{Missing,Tuple{Float64,Float64}},test::Union{Missing,Int64},test2::Union{Missing,Int64}}
        @test all(extract(st, [missing, (9.0, 0.1), (10.0, 0.2), (10.0, 0.3)]) .=== T[
            (geometry = missing, test = missing, test2 = missing)
            (geometry = (9.0, 0.1), test = 1, test2 = 5)
            (geometry = (10.0, 0.2), test = 4, test2 = 8)
            (geometry = (10.0, 0.3), test = missing, test2 = missing)
        ])
        @test extract(st, [missing, (9.0, 0.1), (10.0, 0.2), (10.0, 0.3)]; skipmissing=true) == [
            (geometry = (10.0, 0.2), test = 4, test2 = 8)
        ]
        @test extract(st2, [missing, (2, 2), (2, 1)]; skipmissing=true) == [
            (geometry = (2, 1), a = 7.0, b = 2.0)
        ]
        @test extract(st, [missing, (9.0, 0.1), (10.0, 0.2), (10.0, 0.3)]; skipmissing=true, geometry=false) == [
            (test = 4, test2 = 8)
        ]
        T = @NamedTuple{index::Union{Missing, Tuple{Int,Int}}, test::Union{Missing, Int64}, test2::Union{Missing, Int64}}
        @test extract(st, [missing, (9.0, 0.1), (10.0, 0.2), (10.0, 0.3)]; skipmissing=true, geometry=false, index=true) == T[
            (index = (2, 2), test = 4, test2 = 8)
        ]
        # Subset with `names`
        T = @NamedTuple{geometry::Union{Missing, Tuple{Float64, Float64}}, test2::Union{Missing, Int64}}
        @test all(extract(st, [missing, (9.0, 0.1), (10.0, 0.2), (10.0, 0.3)]; name=(:test2,)) .=== T[
            (geometry = missing, test2 = missing)
            (geometry = (9.0, 0.1), test2 = 5)
            (geometry = (10.0, 0.2), test2 = 8)
            (geometry = (10.0, 0.3), test2 = missing)
        ])
        # Subset with `names` and `skipmissing` with mixed missingvals
        @test extract(st, [missing, (9.0, 0.1), (10.0, 0.2), (10.0, 0.3)]; name=(:test2,), skipmissing=true) == [
            (geometry = (10.0, 0.2), test2 = 8)
        ]
        @test extract(st, [missing, (9.0, 0.1), (10.0, 0.2), (10.0, 0.3)]; name=(:test,), skipmissing=true) == [
            (geometry = (9.0, 0.1), test = 1)
            (geometry = (10.0, 0.2), test = 4)
        ]
    end

    @testset "Errors" begin
        # Missing coord errors
        @test_throws ArgumentError extract(rast, [(0.0, missing), (9.0, 0.1), (9.0, 0.2), (10.0, 0.3)])
        @test_throws ArgumentError extract(rast, [(9.0, 0.1), (0.0, missing), (9.0, 0.2), (10.0, 0.3)])
        @test_throws ArgumentError extract(rast, [(X=0.0, Y=missing), (9.0, 0.1), (9.0, 0.2), (10.0, 0.3)])
    end
end

@testset "Polygons" begin
    # Extract a polygon
    T = @NamedTuple{geometry::Union{Tuple{Float64,Float64}},test::Union{Missing,Int64}}
    @test all(extract(rast_m, poly) .=== T[
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (10.0, 0.1), test = 3)
        (geometry = (10.0, 0.2), test = missing)
    ])
    # Test all the keyword combinations
    @test all(extract(rast_m, poly) .=== T[
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (10.0, 0.1), test = 3)
        (geometry = (10.0, 0.2), test = missing)
    ])
    T = @NamedTuple{test::Union{Missing,Int64}}
    @test all(extract(rast_m, poly; geometry=false) .=== T[
        (test = 1,)
        (test = 3,)
        (test = missing,)
    ])
    T = @NamedTuple{index::Union{Missing,Tuple{Int,Int}},test::Union{Missing,Int64}}
    @test all(extract(rast_m, poly; geometry=false, index=true) .=== T[
        (index = (1, 1), test = 1)
        (index = (2, 1), test = 3)
        (index = (2, 2), test = missing)
    ])
    T = @NamedTuple{geometry::Union{Tuple{Float64,Float64}},index::Union{Missing,Tuple{Int,Int}},test::Union{Missing,Int64}}
    @test all(extract(rast_m, poly; index=true) .=== T[
         (geometry = (9.0, 0.1), index = (1, 1), test = 1)
         (geometry = (10.0, 0.1), index = (2, 1), test = 3)
         (geometry = (10.0, 0.2), index = (2, 2), test = missing)
    ])
    @test extract(rast_m, poly; skipmissing=true) == [
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (10.0, 0.1), test = 3)
    ]                                                         
    @test extract(rast_m, poly; skipmissing=true, geometry=false) == [
        (test = 1,)
        (test = 3,)
    ]                                                         
    @test extract(rast_m, poly; skipmissing=true, geometry=false, index=true) == [
        (index = (1, 1), test = 1)
        (index = (2, 1), test = 3)
    ]                                                         
    @test extract(rast_m, poly; skipmissing=true, index=true) == [
        (geometry = (9.0, 0.1), index = (1, 1), test = 1)
        (geometry = (10.0, 0.1), index = (2, 1), test = 3)
    ]          
    @test extract(rast2, poly; skipmissing=true) == [
        (geometry = (10.0, 0.1), test2 = 7)
        (geometry = (10.0, 0.2), test2 = 8)
    ]                                               
    T = @NamedTuple{geometry::Union{Tuple{Float64,Float64}},test::Union{Missing,Int64}}
    @test all(extract(rast_m, poly) .=== T[
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (10.0, 0.1), test = 3)
        (geometry = (10.0, 0.2), test = missing)
    ])

    @testset "Vector of polygons" begin
        ex = extract(rast_m, [poly, poly, poly])
        @test eltype(ex) == T
        @test all(ex .=== T[
            (geometry = (9.0, 0.1), test = 1)
            (geometry = (10.0, 0.1), test = 3)
            (geometry = (10.0, 0.2), test = missing)
            (geometry = (9.0, 0.1), test = 1)
            (geometry = (10.0, 0.1), test = 3)
            (geometry = (10.0, 0.2), test = missing)
        ])
    end

    @testset "LineString" begin
        @test extract(rast, linestring) == T[
            (geometry=(9.0, 0.1), test=1)
            (geometry=(10.0, 0.1), test=3)
        ]
    end
end

@testset "Extract a linestring" begin
    T = @NamedTuple{geometry::Tuple{Float64,Float64},test::Union{Missing,Int64}}
    Tsm = @NamedTuple{geometry::Tuple{Float64,Float64},test::Int64}

    @test all(extract(rast, l) .=== T[
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (10.0, 0.1), test = 3)
        (geometry = (10.0, 0.2), test = missing)
    ])
end

@testset "Table" begin
    T = @NamedTuple{geometry::Union{Missing, Tuple{Float64, Float64}}, test::Union{Missing, Int64}}
    @test all(extract(rast, table) .=== T[
        (geometry = missing, test = missing)
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (9.0, 0.2), test = 2)
        (geometry = (10.0, 0.3), test = missing)
        (geometry = (10.0, 0.2), test = 4)
    ])
    @test extract(rast, table; skipmissing=true) == [
        (geometry = (9.0, 0.1), test = 1)
        (geometry = (9.0, 0.2), test = 2)
        (geometry = (10.0, 0.2), test = 4)
    ]
    @test extract(rast, table; skipmissing=true, geometry=false) == [
        (test = 1,)
        (test = 2,)
        (test = 4,)
    ]
    @test extract(rast, table; skipmissing=true, geometry=false, index=true) == [
        (index = (1, 1), test = 1,)
        (index = (1, 2), test = 2,)
        (index = (2, 2), test = 4,)
    ]
    @test_throws ArgumentError extract(rast, (foo = zeros(4),))
end

@testset "Empty geoms" begin
    @test extract(rast, []) == NamedTuple{(:geometry, :test),Tuple{Missing,Missing}}[]
    @test extract(rast, []; geometry=false) == NamedTuple{(:test,),Tuple{Missing}}[]
end

#= =#
using BenchmarkTools
using ProfileView

dimz = (X(5.0:0.1:15.0; sampling=Intervals(Start())), Y(-0.1:0.01:0.5; sampling=Intervals(Start())))
rast_m = Raster(rand([missing, 1, 2], dimz); name=:test, missingval=missing)
rast = Raster(rand(Int, dimz); name=:test2, missingval=5)
st = RasterStack((a=rast, b=rast, c=Float32.(rast)))

ks = ntuple(x -> gensym(), 100)
layers = NamedTuple{ks}(ntuple(x -> rast2, length(ks)))
st100 = RasterStack(layers)

poly = GI.Polygon([[(8.0, 0.0), (11.0, 0.0), (11.0, 0.4), (8.0, 0.0)]])
linestring = GI.LineString([(8.0, 0.0), (11.0, 0.0), (11.0, 0.4), (8.0, 0.0)])
line = GI.Line([(8.0, 0.0), (12.0, 4.0)])
polys = [poly for i in 1:2]
lines = [line for i in 1:2]
linestrings = [linestring for i in 1:2]

extract(rast, linestring)
extract(rast, poly)
@time extract(rast, lines; geometry=false, threaded=true, flatten=true)
@time extract(st, lines; geometry=false, threaded=true, flatten=false)
@time extract(st, lines; geometry=false, threaded=false, flatten=false) |> length
@time extract(rast_m, lines; geometry=false, threaded=true, flatten=false) |> length
@time extract(st, lines; geometry=false, threaded=false, flatten=false) |> length
@time extract(st100, lines; geometry=false, threaded=false, flatten=false) |> length
@time extract(st100, lines; geometry=false, threaded=true, flatten=false) |> length
@time extract(rast, lines; geometry=false, threaded=true, flatten=true)
@time extract(rast, lines; threaded=true, flatten=true, progress=false)
@time extract(rast, lines; threaded=false, flatten=false, progress=false);
@benchmark extract(rast, lines; threaded=false, flatten=true)
@benchmark extract(rast, lines; threaded=true, flatten=true)
@benchmark extract(rast, lines; threaded=false, flatten=false)
@benchmark extract(rast, lines; threaded=true, flatten=false)
@profview 
@time extract(rast, polys; flatten=false);
extract(rast, polys; flatten=true);
extract(st100, linestring)
extract(st100, poly)

@b extract(rast2, linestring)
@b extract(rast, linestring)
@b extract(rast, poly)
@b extract(rast, polies)
@b extract(rast, polies; flatten=false)
@b extract(st26, linestring)
@b extract(st26, poly)


f(rast, ls, n; skipmissing=true) = for _ in 1:n extract(rast, ls; skipmissing) end
@profview f(rast2, polies, 10000; skipmissing=false)
@profview f(rast2, poly, 10000; skipmissing=false)
@profview f(rast, poly, 10000; skipmissing=false)
@profview f(rast2, linestring, 10000; skipmissing=false)
@profview f(rast, linestring, 10000; skipmissing=false)
@profview f(rast2, linestring, 10000)
@profview f(st, linestring, 100000; skipmissing=false)
@profview f(st, linestring, 100000; skipmissing=true)
@profview f(st2, linestring, 100000)
@profview f(st26, linestring, 1000)
@profview f(st26, line, 1000)
@profview f(st26, poly, 10000; skipmissing=true)
@profview f(st26, linestring, 1000; skipmissing=true)
@profview f(st26, linestring, 1000; skipmissing=false)

# Demo plots
using GLMakie

Makie.plot(rast; colormap=:Reds)
Makie.plot!(polys)
points = getindex.(extract(rast, poly; index=true, flatten=true), :geometry)
Makie.scatter!(points; color=(:yellow, 0.5))

Makie.plot(rast; colormap=:Reds)
Makie.plot!(polys)
points = getindex.(extract(rast, poly; index=true, flatten=true, boundary=:touches), :geometry)
Makie.scatter!(points; color=(:yellow, 0.5))

Makie.plot(rast; colormap=:Reds)
Makie.plot!(linestring; color=:violet, linewidth=5)
points = getindex.(extract(rast, linestring; index=true), :geometry)
Makie.scatter!(points; color=:yellow)


#= =#
