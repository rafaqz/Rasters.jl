using Rasters, Extents, ArchGDAL, Test
import GeoInterface as GI

a = Raster((1:26) * (1:31)', (X(-20:5), Y(0:30)))
b = a .* 2
c = cat(a .* 3, a .* 3; dims=:newdim)
st = RasterStack((; a, b, c))
a_m = let 
    a = Array{Union{Missing,Int}}(undef, 26, 31)
    a .= (1:26) * (1:31)'
    a[1:10, 3:10] .= missing
    Raster(a, (X(-20:5), Y(0:30)))
end

pointvec = [(-20.0, 30.0), (-20.0, 10.0), (0.0, 10.0), (0.0, 30.0), (-20.0, 30.0)]
pointvec_empty = [(-100.0, 0.0), (-100.0, 0.0), (-100.0, 0.0), (-100.0, 0.0), (-100.0, 0.0)]
pointvec_out_of_bounds = [(-40.0, -40.0), (-40.0, -35.0), (-35.0, -35.0), (-35.0, -40.0)]
linestring = GI.LineString(pointvec)
polygon = GI.Polygon([pointvec])
polygon_empty = GI.Polygon([pointvec_empty])
polygon_out_of_bounds = GI.Polygon([pointvec_out_of_bounds])

@testset "polygons" begin
    @test zonal(sum, a; of=polygon) ==
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

    @test zonal(sum, st; of=polygon) == zonal(sum, st; of=[polygon])[1] ==
        zonal(sum, st; of=(geometry=polygon, x=:a, y=:b)) ==
        zonal(sum, st; of=[(geometry=polygon, x=:a, y=:b)])[1] ==
        zonal(sum, st; of=[(geometry=polygon, x=:a, y=:b)])[1] ==
        maplayers(sum ∘ skipmissing, mask(st; with=polygon))  
    @test zonal(sum, st; of=st) == 
        zonal(sum, st; of=dims(st)) == 
        zonal(sum, st; of=Extents.extent(st)) == 
        sum(st)

    @testset "skipmissing" begin
        @test zonal(sum, a_m; of=polygon, skipmissing=false) === missing
        @test zonal(sum, a_m; of=polygon, skipmissing=true) isa Int
        @test !zonal(x -> x isa Raster, a_m; of=polygon, skipmissing=true)
        @test zonal(x -> x isa Raster, a_m; of=polygon, skipmissing=false)
    end
end


@testset "Lines" begin
    @test 
end

using BenchmarkTools
@btime zonal(sum, a; of=linestring, skipmissing=true)
@btime zonal(x -> sum(x), a; of=linestring, skipmissing=true) 
@btime zonal(sum ∘ skipmissing, a; of=linestring, skipmissing=false)
using ProfileView
f(a, n) = for _ in 1:n zonal(sum ∘ skipmissing, a; of=linestring, skipmissing=false) end
g(a, n) = for _ in 1:n zonal(sum ∘ skipmissing, a; of=linestring, skipmissing=true) end
using Cthulhu
Cthulhu.@descend zonal(sum, a; of=linestring, skipmissing=true)
ProfileView.descend_clicked()
VSCodeServer.@profview 
f(a, 1000)
ProfileView.@profview f(a, 100000)
ProfileView.@profview g(a, 100000)
zonal(a; of=[linestring, linestring], skipmissing=true) do A
    @show A
    sum(A)
end
zonal(a; of=[linestring, linestring], skipmissing=false) do A
    @show collect(skipmissing(A))
    sum(skipmissing(A))
end zonal(sum, a; of=linestring, skipmissing=true)

@testset "return missing" begin
    @test ismissing(zonal(sum, a; of=[polygon, polygon_out_of_bounds, polygon])[2]) &&
        ismissing(zonal(sum, a; of=[polygon_out_of_bounds, polygon])[1]) &&
        ismissing(zonal(sum, a; of=(geometry=polygon_out_of_bounds, x=:a, y=:b))) &&
        ismissing(zonal(sum, a; of=[(geometry=polygon_out_of_bounds, x=:a, y=:b)])[1])
    @test zonal(sum, a; of=[polygon_out_of_bounds, polygon_out_of_bounds, polygon])[3] == 
        sum(skipmissing(mask(a; with=polygon)))
end