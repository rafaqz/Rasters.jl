using Rasters, ArchGDAL
using Test

include(joinpath(dirname(pathof(Rasters)), "../test/test_utils.jl"))
url = "https://download.osgeo.org/geotiff/samples/gdal_eg/cea.tif"
gdalpath = maybedownload(url)

@testset "warp" begin
    # test that warp actually does *something*
    r = Raster(gdalpath)
    crs_ = crs(r).val
    warped = warp(r, Dict(:t_srs => "EPSG:25832"); missingval=nothing)
    @test warped isa Raster
    @test size(warped) == (720, 721)
    # the crs is way off, the image is rotated - all four corners should be black
    missingval(warped) === nothing
    @test warped[1, 1] === warped[1, end] === warped[end, 1] === warped[end, end] === 0x00

    warped = warp(r, Dict(:t_srs => "EPSG:25832"))
    @test warped isa Raster
    @test size(warped) == (720, 721)
    # the crs is way off, the image is rotated - all four corners should be black
    missingval(warped) === nothing
    @test warped[1, 1] === warped[1, end] === warped[end, 1] === warped[end, end] === 0xff
    # now compute mean squared error of the back transformation
    warped_back = Rasters.trim(warp(warped, Dict(:t_srs => crs_), res=map(step, lookup(r))))
    # subtracting UInts brings us into hell -> Int
    # we also need to shrink the range because of some bleed during warp
    diff_ = Int.(warped_back[2:end-1, 2:end-1]) .- r
    @test sum(x->x^2, diff_) / prod(size(diff_)) < 600
end
