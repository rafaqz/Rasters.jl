using Test
using Rasters
using DimensionalData
using DimensionalData.Lookups
using DimensionalData.Dimensions
using Dates
using CommonDataModel
using ZarrDatasets
using Zarr
using Rasters: FileArray, FileStack, Zarrsource, crs, bounds, name, trim

# A small local store shaped like a 2.5 degree, 16 day data cube, with consolidated metadata
path = tempname() * ".zarr"
cubedims = (X(-178.75:2.5:178.75), Y(-88.75:2.5:88.75), Ti(DateTime(1979, 1, 9):Day(16):DateTime(1979, 6, 18)))
cube(; kw...) = Raster(Array{Union{Missing,Float32}}(rand(Float32, length.(cubedims))), cubedims; missingval=missing, kw...)
t2m = cube(; metadata=Rasters.Metadata{Zarrsource}(Dict("original_name" => "t2m")))
t2m[1, 1, 1] = -28.866226f0
t2m[end, end, end] = missing
write(path, RasterStack((; air_temperature_2m=t2m, snow_sublimation=cube())); missingval=-9999f0)
Zarr.consolidate_metadata(Zarr.DirectoryStore(path), "")
constructor_kw = (; consolidated=true)

zraster = Raster(path; name="air_temperature_2m", constructor_kw)
lazyarray = Raster(path; lazy=true, name="air_temperature_2m", constructor_kw)
eagerarray = Raster(path; lazy=false, name="air_temperature_2m", constructor_kw)
@test_throws ArgumentError Raster("notafile.zarr/")

@testset "Raster from dataset" begin
    ds = ZarrDatasets.ZarrDataset(path; constructor_kw...)
    var = CommonDataModel.variable(ds, "air_temperature_2m")
    dsarray = Raster(ds; name=:air_temperature_2m)
    vararray = Raster(var; name=:air_temperature_2m)
    @test dims(dsarray) == dims(vararray) == dims(zraster)
    @test size(dsarray) == size(vararray) == size(zraster)
    @test all(dsarray .=== vararray  .=== zraster)
end

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
    @time read!(zraster, A2);
    A3 = copy(A) .= 0
    @time read!(zraster, A3)
    @test all(A .=== A2) 
    @test all(A .=== A3)
end

@testset "array properties" begin
    @test name.(dims(zraster)) == (:X, :Y, :Ti)
    @test length(dims(zraster, X)) == 144
    @test lookup(zraster,X) == collect(-178.75:2.5:178.75)
    # TODO the spatial bounds are strange, because the data is point data
    # We should find a dataset that has actual intervals

    @test bounds(zraster) == (
        (-178.75, 178.75),
        (-88.75, 88.75),
        (DateTime("1979-01-09T00:00:00"), DateTime("1979-06-18T00:00:00")),
    )
end

@testset "dimensions" begin
    @test ndims(zraster) == 3
    @test length.(dims(zraster)) == (144, 72, 11)
    @test dims(zraster) isa Tuple{<:X,<:Y,<:Ti}
    @test refdims(zraster) == ()
    @test val.(span(zraster)) == (2.5, 2.5, (nothing, nothing))
    @test typeof(lookup(zraster)) <: Tuple{<:Mapped,<:Mapped,<:Sampled}
end

@testset "other fields" begin
    @test ismissing(missingval(zraster))
    @test metadata(zraster)["original_name"] == "t2m"
    @test metadata(zraster) isa Rasters.Metadata{<:Rasters.CDMsource, Dict{String, Any}}
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

@testset "CF conventions" begin
    path = tempname() * ".zarr"
    global_attrs = Dict(
        "Conventions" => "CF-1.7",
        "_FillValue" => -999.0,
        "units" => "K",
        "project" => "GOES",
    )
    zg = zgroup(path; attrs = global_attrs)
    x_attrs = Dict(
        "units"             => "rad",
        "add_offset"        => 0,
        "_ARRAY_DIMENSIONS" => ["x"],
        "axis"              => "X",
        "long_name"         => "GOES fixed grid projection x-coordinate",
        "scale_factor"      => 1,
        "_Netcdf4Dimid"     => 1,
        "standard_name"     => "projection_x_coordinate",
        "NAME"              => "x",
    )
    xs = zcreate(Float64, zg, "x", 100; attrs = x_attrs)
    xs .= LinRange(-5.434177815866809e6, 5.322929839864517e6, 100)

    y_attrs = Dict(
        "units"             => "rad",
        "add_offset"        => 5.258454335331338e6,
        "_ARRAY_DIMENSIONS" => ["y"],
        "axis"              => "Y",
        "long_name"         => "GOES fixed grid projection y-coordinate",
        "scale_factor"      => -1.0625617981243342e7,
        "_Netcdf4Dimid"     => 1,
        "standard_name"     => "projection_y_coordinate",
        "NAME"              => "y",
    )
    ys = zcreate(Float64, zg, "y", 100; attrs = y_attrs)
    ys .= LinRange(0, 1, 100)

    val_attrs = Dict(
        "_Netcdf4Dimid" => 0,
        "scale_factor" => 2.0,
        "add_offset" => 180.0,
        "_FillValue" => -1, 
        "units" => "K", 
        "_ARRAY_DIMENSIONS" => ["y", "x"], 
        "grid_mapping" => "goes_imager_projection",
        "valid_range" => [0, -6],
    )
    vals = zcreate(Float64, zg, "values", 100, 100; attrs = val_attrs)
    vals .= (data = rand(100, 100))
    vals[1, 1] = 1.0
    vals[end, end] = 0.0

    ra = Raster(path)

    @test extrema(ra) == (180.0, 182.0) # test scale and offset
    @test Rasters.GeoInterface.extent(ra) ≈ Rasters.GeoInterface.Extent(X = (-5.434177815866809e6, 5.322929839864517e6), Y = (-5.367163645912004e6, 5.258454335331338e6))
    @test Rasters.isreverse(ra.dims[2])
    @test Rasters.isforward(ra.dims[1])
    @test collect(extrema(ra.dims[1])) ≈ collect(extrema(xs))
    @test collect(extrema(ra.dims[2])) ≈ collect(reverse(extrema(ys)) .* y_attrs["scale_factor"] .+ y_attrs["add_offset"])

    @testset "write" begin
        fn = tempname() * ".zarr"
        write(fn, ra)
        @test all(Raster(fn) .=== ra)
        # Broken: Rasters opens for writing with mode "a", which ZarrDatasets rejects
        x = Raster(fn; lazy=true)
        @test_broken (open(O -> O .= 1, x; write=true); all(Raster(fn) .== 1))
    end
end


zarrstack = RasterStack(path; lazy=true, constructor_kw)
@testset "RasterStack" begin
    @test all(zarrstack.snow_sublimation .=== Raster(path; name=:snow_sublimation, constructor_kw))
end
@testset "RasterStack from dataset" begin
    ds = ZarrDatasets.ZarrDataset(path; constructor_kw...)
    dsstack = RasterStack(ds; lazy=true)
    @test dims(dsstack) == dims(zarrstack)
    @test size(dsstack) == size(zarrstack)
end

@testset "backend-specific open kwargs persist across lazy reopen" begin
    # Regression test: a lazily-read layer that is only discoverable under a
    # backend-specific open option (e.g. Zarr `consolidated`) must still be
    # reachable once the underlying `FileArray` reopens the store to read data,
    # not just at initial layer discovery.
    #
    # Simulates a real-world case (ARCO-ERA5's incomplete store listing, fixed
    # by ZarrDatasets.jl's `consolidated` keyword): consolidate metadata while
    # both `a` and `b` exist, then delete `b`'s `.zarray` marker so plain
    # listing can no longer find it -- but its chunk data and consolidated
    # metadata entry remain, so `consolidated=true` can still discover and
    # read it.
    kwpath = tempname() * ".zarr"
    kwzg = zgroup(kwpath)
    ka = zcreate(Float64, kwzg, "a", 4; attrs=Dict("_ARRAY_DIMENSIONS" => ["x"]))
    ka .= [1.0, 2.0, 3.0, 4.0]
    kb = zcreate(Float64, kwzg, "b", 4; attrs=Dict("_ARRAY_DIMENSIONS" => ["x"]))
    kb .= [10.0, 20.0, 30.0, 40.0]
    Zarr.consolidate_metadata(kwzg.storage, "")
    rm(joinpath(kwpath, "b", ".zarray"))

    # Without the open option, `b` isn't visible at all.
    st_plain = RasterStack(kwpath; source=Zarrsource(), lazy=true)
    @test propertynames(st_plain) == (:a,)

    # With it, `b` is discovered at open time...
    st_cons = RasterStack(kwpath; source=Zarrsource(), lazy=true, constructor_kw=(; consolidated=true))
    @test :b in propertynames(st_cons)

    # ...and reading it -- which reopens the store via `FileArray` -- must
    # still find it, i.e. `constructor_kw` has to be reapplied, not just used once.
    @test collect(st_cons[:b]) == [10.0, 20.0, 30.0, 40.0]
    @test collect(st_cons[:a]) == [1.0, 2.0, 3.0, 4.0]

    # Unrecognized top-level keywords must still error -- `constructor_kw` is a
    # deliberately separate channel, not a generic catch-all that would mask
    # keyword typos/misuse.
    @test_throws MethodError RasterStack(kwpath; source=Zarrsource(), lazy=true, totally_bogus_kwarg=true)
end
