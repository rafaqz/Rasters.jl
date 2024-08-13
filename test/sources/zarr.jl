using Rasters
using Zarr
using ZarrDatasets
using Rasters: FileArray, FileStack, Zarrsource, crs, bounds, name, trim

path = "https://s3.bgc-jena.mpg.de:9000/esdl-esdc-v3.0.2/esdc-16d-2.5deg-46x72x1440-3.0.2.zarr"

@testset "Zarr Raster open" begin

zraster = Raster(path; name="air_temperature_2m")
lazyarray = Raster(path; lazy=true, name="air_temperature_2m")
eagerarray = Raster(path; lazy=false, name="air_temperature_2m")
@test_throws ArgumentError Raster("notafile.zarr/")

@testset "lazyness" begin
    # Eager is the default
    @test parent(zraster) isa Array
    @test parent(lazyarray) isa FileArray
    @test parent(eagerarray) isa Array
end
@testset "read" begin
    @time A = read(lazyarray);
    @test A isa Raster
    @test parent(A) isa Array
    A2 = copy(A) .= 0
    @time read!(ncarray, A2);
    A3 = copy(A) .= 0
    @time read!(ncsingle, A3)
    @test all(A .=== A2) 
    @test all(A .=== A3)
end

@testset "array properties" begin
    @test name.(dims(zraster)) == (:X, :Y, :Ti)
    @test length(dims(zraster, X)) == 144
    @test index(zraster,X) == collect(-178.75:2.5:178.75)
    # TODO the spatial bounds are strange, because the data is point data
    # We should find a dataset that has actual intervals
    @test bounds(zraster) == (
        (-178.75, 178.75),
        (-88.75, 88,75),
        (DateTime("1979-01-09T00:00:00"), DateTime("2021-12-27T00:00:00")),
    )
end
@testset "dimensions" begin
    @test ndims(zraster) == 3
    @test length.(dims(zraster)) == (144, 72, 989)
    @test dims(zraster) isa Tuple{<:X,<:Y,<:Ti}
    @test refdims(zraster) == ()
    @test val.(span(ncarray)) == (2.5, 2.5, (nothing, nothing))
    @test typeof(lookup(ncarray)) <: Tuple{<:Mapped,<:Mapped,<:Sampled}
end
@testset "other fields" begin
    @test ismissing(missingval(zraster))
    @test metadata(r)["original_name"] == "t2m"
    @test metadata(zraster) isa Metadata{<:Rasters.CDMsource, Dict{String, Any}}
    @test name(zraster) == :air_temperature_2m
end

@testset "indexing" begin
    @test zraster[Ti(1)] isa Raster{<:Any,2}
    @test zraster[Y(1), Ti(1)] isa Raster{<:Any,1}
    @test zraster[X(1), Ti(1)] isa Raster{<:Any,1}
    @test zraster[X(1), Y(1), Ti(1)] == -28.866226f0 == parent(zraster)[1,1,1]
    @test zraster[X(30), Y(30), Ti(1)] isa Float32
    # Alaska
    @test zraster[Y(Near(-88.75)), X(Near(-178.74)), Ti(1)] ==-28.866226f0
    @test zraster[Ti(At(DateTime(1979,1,9))), X(At(-178.75)), Y(At(-88.75))] == -28.866226f0
end

end
