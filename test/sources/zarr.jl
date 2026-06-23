using Test
using Rasters
using DimensionalData
using DimensionalData.Lookups
using DimensionalData.Dimensions
using Dates
using CommonDataModel
using ZarrDatasets
using ZarrDatasets.Zarr
using Rasters: FileArray, FileStack, Zarrsource, crs, bounds, name, trim

# Build a small CF-compliant Zarr group on disk so these tests run offline,
# without depending on a remote dataset being reachable.
const ZARR_NX, ZARR_NY, ZARR_NT = 10, 8, 5
const ZARR_LON = collect(-178.75 .+ 2.5 .* (0:ZARR_NX-1))
const ZARR_LAT = collect(-88.75 .+ 2.5 .* (0:ZARR_NY-1))
# Deterministic data so cell values are predictable: value = x + 10y + 100t
const ZARR_DATA = Float32[i + 10j + 100k for i in 1:ZARR_NX, j in 1:ZARR_NY, k in 1:ZARR_NT]

function build_test_zarr(path)
    zg = zgroup(path; attrs=Dict("Conventions" => "CF-1.7"))
    lon = zcreate(Float64, zg, "lon", ZARR_NX;
        attrs=Dict("_ARRAY_DIMENSIONS" => ["lon"], "axis" => "X", "standard_name" => "longitude"))
    lon[:] = ZARR_LON
    lat = zcreate(Float64, zg, "lat", ZARR_NY;
        attrs=Dict("_ARRAY_DIMENSIONS" => ["lat"], "axis" => "Y", "standard_name" => "latitude"))
    lat[:] = ZARR_LAT
    time = zcreate(Float64, zg, "time", ZARR_NT;
        attrs=Dict("_ARRAY_DIMENSIONS" => ["time"], "axis" => "T", "standard_name" => "time",
            "units" => "days since 1979-01-09", "calendar" => "standard"))
    time[:] = collect(16.0 .* (0:ZARR_NT-1))
    # Zarr stores `_ARRAY_DIMENSIONS` in reverse of the Julia array order
    air = zcreate(Float32, zg, "air_temperature_2m", ZARR_NX, ZARR_NY, ZARR_NT;
        attrs=Dict("_ARRAY_DIMENSIONS" => ["time", "lat", "lon"],
            "_FillValue" => -9999.0f0, "original_name" => "t2m"))
    air[:, :, :] = ZARR_DATA
    snow = zcreate(Float32, zg, "snow_sublimation", ZARR_NX, ZARR_NY, ZARR_NT;
        attrs=Dict("_ARRAY_DIMENSIONS" => ["time", "lat", "lon"]))
    snow[:, :, :] = ZARR_DATA .* 2
    return path
end

path = tempname() * ".zarr"
build_test_zarr(path)

zraster = Raster(path; name="air_temperature_2m")
lazyarray = Raster(path; lazy=true, name="air_temperature_2m")
eagerarray = Raster(path; lazy=false, name="air_temperature_2m")
@test_throws ArgumentError Raster("notafile.zarr/")

@testset "Raster from dataset" begin
    ds = ZarrDatasets.ZarrDataset(path)
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
    A = read(lazyarray);
    @test A isa Raster
    @test parent(A) isa Array
    A2 = copy(A) .= 0
    read!(zraster, A2);
    A3 = copy(A) .= 0
    read!(zraster, A3)
    @test all(A .=== A2)
    @test all(A .=== A3)
end

@testset "array properties" begin
    @test name.(dims(zraster)) == (:X, :Y, :Ti)
    @test length(dims(zraster, X)) == ZARR_NX
    @test lookup(zraster, X) == ZARR_LON
    @test lookup(zraster, Y) == ZARR_LAT
end

@testset "dimensions" begin
    @test ndims(zraster) == 3
    @test length.(dims(zraster)) == (ZARR_NX, ZARR_NY, ZARR_NT)
    @test dims(zraster) isa Tuple{<:X,<:Y,<:Ti}
    @test refdims(zraster) == ()
    @test bounds(zraster) == (
        (-178.75, -156.25),
        (-88.75, -71.25),
        (DateTime("1979-01-09T00:00:00"), DateTime("1979-03-14T00:00:00")),
    )
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
    # Data is value = x + 10y + 100t for 1-based indices
    @test zraster[X(1), Y(1), Ti(1)] == 111f0 == parent(zraster)[1, 1, 1]
    @test zraster[X(2), Y(1), Ti(1)] == 112f0
    @test zraster[X(1), Y(2), Ti(1)] == 121f0
    @test zraster[X(1), Y(1), Ti(2)] == 211f0
    @test zraster[Y(Near(-88.75)), X(Near(-178.74)), Ti(1)] == 111f0
    @test zraster[Ti(At(DateTime(1979, 1, 9))), X(At(-178.75)), Y(At(-88.75))] == 111f0
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

    ext = Rasters.GeoInterface.extent(ra)
    @test extrema(ra) == (180.0, 182.0) # test scale and offset
    @test all(ext.X .≈ (-5.434177815866809e6, 5.322929839864517e6))
    @test all(ext.Y .≈ (-5.367163645912004e6, 5.258454335331338e6))
    @test Rasters.isreverse(ra.dims[2])
    @test Rasters.isforward(ra.dims[1])
    @test all(extrema(ra.dims[1]) .≈ extrema(xs))
    @test all(extrema(ra.dims[2]) .≈ reverse(extrema(ys)) .* y_attrs["scale_factor"] .+ y_attrs["add_offset"])

    @testset "write" begin
        fn = tempname() * ".zarr"
        write(fn, ra)
        @test all(Raster(fn) .=== ra)
        # In-place writing to an existing Zarr needs append mode, which
        # ZarrDatasets.jl does not yet support.
        x = Raster(fn; lazy=true)
        @test_broken begin
            open(x; write=true) do O
                O .= 1
            end
            all(Raster(fn) .== 1)
        end
    end

    @testset "write with force" begin
        # Issue #1060: force=true should overwrite an existing zarr
        fn = tempname() * ".zarr"
        write(fn, ra)
        @test isdir(fn)
        # Without force should error
        @test_throws ArgumentError write(fn, ra)
        # With force should overwrite
        write(fn, ra; force=true)
        @test all(Raster(fn) .=== ra)
        # Overwrite with different data
        ra2 = rebuild(ra; data=ones(size(ra)))
        write(fn, ra2; force=true)
        @test all(Raster(fn) .=== ra2)
    end
end


zarrstack = RasterStack(path; lazy=true)
@testset "RasterStack" begin
    @test all(zarrstack.snow_sublimation .=== Raster(path; name=:snow_sublimation))
end
@testset "RasterStack from dataset" begin
    ds = ZarrDatasets.ZarrDataset(path)
    dsstack = RasterStack(ds; lazy=true)
    @test dims(dsstack) == dims(zarrstack)
    @test size(dsstack) == size(zarrstack)
end
