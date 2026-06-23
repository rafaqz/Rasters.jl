using Rasters, Test, Statistics
using Rasters: FileArray, FileStack
using Rasters.DimensionalData: X, Y
using Rasters.DimensionalData
using NCDatasets
using ArchGDAL

const ext = ".tif"

@testset "open/close — GDAL Raster" begin
    A = Raster(rand(Float32, 5, 5), (X(1:5), Y(1:5)); name=:test)
    filename = tempname() * ext
    write(filename, A; force=true)
    r = Raster(filename; lazy=true)

    @test parent(r) isa AbstractArray

    o = open(r)
    @test parent(o) !== parent(r)            # parent was swapped for opened var
    @test o[1:1, 1:1] == r[1:1, 1:1]
    @test sum(o) ≈ open(sum, r)              # matches closure form

    close(o)
    @test_throws Exception o[1:1, 1:1]       # dataset released

    # Closure form still works after explicit close
    @test open(sum, r) isa Float32
end

@testset "open/close — GDAL write" begin
    A = Raster(rand(Float32, 5, 5), (X(1:5), Y(1:5)); name=:test)
    filename = tempname() * ext
    write(filename, A; force=true)
    r = Raster(filename; lazy=true)

    before = read(r)
    o = open(r; write=true)
    o .*= 2
    close(o)
    @test read(Raster(filename)) ≈ before .* 2
end

@testset "open is a no-op for in-memory rasters" begin
    r = Raster(rand(Float32, 3, 3), (X(1:3), Y(1:3)))
    @test open(r) === r
    @test close(r) === r
end

@testset "open/close — NetCDF Raster" begin
    filename = tempname() * ".nc"
    A = Raster(rand(Float32, 4, 4), (X(1:4), Y(1:4)); name=:a)
    write(filename, A; force=true)
    r = Raster(filename; lazy=true)

    o = open(r)
    @test o[1:1, 1:1] == r[1:1, 1:1]
    @test sum(o) ≈ open(sum, r)

    close(o)
    @test_throws Exception o[1:1, 1:1]
end

@testset "open/close — NetCDF Stack" begin
    filename = tempname() * ".nc"
    A1 = Raster(rand(Float32, 4, 4), (X(1:4), Y(1:4)); name=:a)
    A2 = Raster(rand(Float32, 4, 4), (X(1:4), Y(1:4)); name=:b)
    st_mem = RasterStack((a=A1, b=A2))
    write(filename, st_mem; force=true)
    st = RasterStack(filename; lazy=true)
    @test parent(st) isa FileStack

    o = open(st)
    @test parent(o) isa Rasters.OpenStack
    @test o[:a][1:1, 1:1] == st[:a][1:1, 1:1]
    @test sum(o[:a]) ≈ open(s -> sum(s[:a]), st)

    close(o)
    @test_throws Exception o[:a][1:1, 1:1]

    # Closure form still works
    @test open(s -> sum(s[:a]), st) isa Union{Missing,Float32}
end
