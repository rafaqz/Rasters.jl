using Rasters, ArchGDAL
using Test

include(joinpath(dirname(pathof(Rasters)), "../test/test_utils.jl"))
url = "https://download.osgeo.org/geotiff/samples/gdal_eg/cea.tif"
gdalpath = maybedownload(url)

@testset "warp" begin
    # test that warp actually does *something*
    r = Raster(gdalpath)
    crs_ = crs(r).val
    warped = warp(r, Dict(:t_srs => "EPSG:25832"); missingval=0xff)
    @test warped isa Raster
    @test size(warped) == (720, 721)
    # the crs is rotated so the image is rotated an all four corners should be missing
    missingval(warped) === 0xff
    parent(warped)
    @test warped[1, 1] === warped[1, end] === warped[end, 1] === warped[end, end] === 0xff == missingval(warped)
    # now compute mean squared error of the back transformation
    res = map(step, lookup(r))
    warped_back = Rasters.trim(warp(warped, Dict(:t_srs => crs_); res, missingval=0xff))
    # subtracting UInts brings us into hell -> Int
    # we also need to shrink the range because of some bleed during warp
    diff_ = parent(warped_back[2:end-1, 2:end-1]) .- r

    @test sum(x -> x^2, diff_) / prod(size(diff_)) < 600
end
