using Rasters
import NCDatasets
import NCDatasets.NetCDF_jll
using NearestNeighbors
using OrderedCollections
using Rasters.Lookups
using DimensionalData: MergedLookup
using Test
using Rasters: name, bounds
using Dates
using GeoInterface

#
#using NCDatasets
# using NCDatasets.NetCDF_jll

# NetCDF_jll.ncdump() do exe
#     run(`$exe clim.nc`)
# end
# ds = collect(NCDataset("clim.nc")["time"])
# ds = collect(NCDataset("clim.nc")["climatology_bounds"])

# Build all ncgen files
testdir = joinpath(dirname(dirname(Base.pathof(Rasters))), "test")
cfdir = joinpath(testdir, "cf")
cfdatadir = joinpath(testdir, "data", "cf")
mkpath(cfdatadir)
example_paths = map(filter(endswith(".ncgen"), readdir(cfdir))) do input_name
    input_path = realpath(joinpath(cfdir, input_name))
    output_path = joinpath(cfdatadir, splitext(input_name)[1] * ".nc")
    NetCDF_jll.ncgen() do exe
        run(`$exe -k nc4 -b -o $output_path $input_path`)
    end
    output_path
end

example_ids = replace.(Base.strip.(first.(split.(basename.(example_paths), r"[A-Za-z]")), '_'), ("_" => ".",))
examples = OrderedDict(example_ids .=> example_paths)

rasterstacks = map(enumerate(example_paths)) do (i, example_path)
    name = splitext(basename(example_path))[1]
    println(i, " => ", name)
    name => RasterStack(example_path; lazy=true)
end |> OrderedDict

# rasters = map(example_paths) do example_path
#     name = splitext(basename(example_path))[1]
#     name == "5_1_independent_coordinate_variables" && return name => nothing
#     name => Raster(example_path; lazy=true)
# end
# st = last.(rasterstacks)[12]
# st = last.(rasterstacks)[14]
# refdims(st)
# inds = findall(.!(isempty.(refdims.(last.(rasterstacks)))))
# refdims(last.(rasterstacks[inds])[1])

@testset "2.1 string variable representation" begin
    rast = RasterStack(examples["2.1"]; lazy=true)
    @test rast.char_variable == rast.str_variable
    @test RasterStack(examples["2.1"]; lazy=true).char_variable == rast.str_variable
end

@testset "3.1 units metadata" begin
    rast = RasterStack(examples["3.1"])
    md_os = metadata(rast[:Tonscale])
    md_diff = metadata(rast[:Tdifference])
    @test md_os["units"] == md_diff["units"] == "degC"
    @test md_os["cell_methods"] == md_diff["cell_methods"] == "area: mean"
    @test md_os["standard_name"] == md_diff["standard_name"] == "surface_temperature"

    @test md_os["long_name"] == "global-mean surface temperature"
    @test md_diff["long_name"] == "change in global-mean surface temperature relative to pre-industrial"
    @test md_os["units_metadata"] == "temperature: on_scale"
    @test md_diff["units_metadata"] == "temperature: difference"
end

@testset "5.1 independent coordinate variables" begin
    rast = RasterStack(examples["5.1"]);
    testdims = (X(Sampled(1.0f0:36.0f0)), 
                Y(Sampled(1.0f0:18.0f0)), 
                Dim{:pres}(Sampled(1.0f0:15.0f0)), 
                Ti(Sampled(collect(DateTime(1990):Day(1):DateTime(1990, 1, 4)))))
    @test all(map(==, dims(rast), map(DimensionalData.format, testdims)))
    @test metadata(rast.xwind)["long_name"] == "zonal wind"
    @test metadata(rast.xwind)["units"] == "m/s"
    @test metadata(dims(rast, :pres))["long_name"] == "pressure"
    @test metadata(dims(rast, :pres))["units"] == "hPa"
    @test metadata(dims(rast, X))["long_name"] == "longitude"
    @test metadata(dims(rast, X))["units"] == "degrees_east"
    @test metadata(dims(rast, Y))["long_name"] == "latitude"
    @test metadata(dims(rast, Y))["units"] == "degrees_north"
    @test metadata(dims(rast, Ti))["long_name"] == "time"
    # This one is kind of redundant, we have converted to DateTime already
    @test metadata(dims(rast, Ti))["units"] == "days since 1990-1-1 0:0:0" ;
end

@testset "5.2 two-dimensional coordinate variables" begin
    rast = RasterStack(examples["5.2"]; lazy=true)
    # TODO probably should have real data for this
    @test rast[X=Near(1.0), Y=Near(2.0), Z=1].T isa Float32
end

@testset "5.3 reduced horizontal grid" begin
    rast = RasterStack(examples["5.3"]; lazy=true)
    @test rast[rgrid=2].PS === 0.2f0
    @test rast[rgrid=At(3.0, 30.0)].PS === 0.3f0
    # Selecting on an inner dim of a MergedLookup with At currently fails:
    # DiskArrays.jl can't getindex with `Base.LogicalIndex{Int64, Vector{Bool}}`.
    @test_broken (rast[X=At(4.0)].PS == Union{Missing,Float32}[0.4f0])
    @test_broken (rast[Y=At(30.0)].PS == Union{Missing,Float32}[0.3f0])
    @test rast[X=At(1.0), Y=At(10.0)] == (PS=0.1f0,)
end

@testset "5.6 rotated pole grid" begin
    rast = RasterStack(examples["5.6"]; lazy=true)
    # Rotated pole dimensions are now typed as X and Y (detected via standard_name grid_longitude/grid_latitude)
    @test lookup(rast, X) isa ProjectedArrayLookup
    @test lookup(rast, Y) isa ProjectedArrayLookup
    @test lookup(rast, Z) == DimensionalData.format(Sampled([100.0f0, 200.0f0]; span=Regular(100.0f0)))
    # Selection by geographic coordinates (lon/lat stored in ProjectedArrayLookup matrix)
    @test rast.temp[X=Near(-0.9), Y=Near(3.0), Z=1] === 3.0f0
    @test rast.temp[X=At(-0.44758424f0), Y=At(2.1773858f0), Z=1] === 2.0f0
end

@testset "5.8-9 WGS 84 latitude_longitude" begin
    # CF latitude_longitude grid mapping is now parsed via CFCoordinateReferenceSystems
    # and returns WellKnownText2 (not EPSG) since the CF files don't contain EPSG codes
    rast = RasterStack(examples["5.8"]; lazy=true)
    @test crs(rast) isa WellKnownText2
    rast = RasterStack(examples["5.9"]; lazy=true)
    @test crs(rast) isa WellKnownText2
end

@testset "5.10 British national grid" begin
    rast = RasterStack(examples["5.10"]; lazy=true);
    @test keys(rast) == (:temp, :pres)
    # Multi-grid_mapping ("crsOSGB: x y crsWGS84: lat lon") is resolved per
    # dim, so crs(stack) returns the primary CRS attached to the spatial
    # dims (crsOSGB, transverse mercator) - matching the single-object
    # expectation of every downstream tool.
    @test crs(rast) isa WellKnownText2
    @test crs(dims(rast, X)) isa WellKnownText2
    @test crs(dims(rast, Y)) isa WellKnownText2
    # The secondary mapping (crsWGS84, scoped to the lat/lon aux coords)
    # becomes mappedcrs on the ProjectedArrayLookup.
    @test lookup(dims(rast, X)).mappedcrs isa WellKnownText2
end

@testset "5.11 WGS 84 WellKnownText2" begin
    rast = RasterStack(examples["5.11"]; lazy=true)
    @test crs(rast) isa WellKnownText2
    # Extent now works correctly
    ext = extent(rast)
    @test ext.X == (0.0, 350.0)
    @test ext.Y == (-85.0, 85.0)
end

@testset "5.14 refdims from scalar coordinates" begin
    rast = RasterStack(examples["5.14"]; lazy=true);
    testdims = (Dim{:atime}([0.0]; span=Regular(0.0)), Dim{:p500}([500.0]; span=Regular(0.0)))
    @test all(map(==, refdims(rast), map(DimensionalData.format, testdims)))
end

@testset "domains" begin
    # TODO: load the dimensions of the domain,
    # even though there are no layers that use them?
    # Its a kind of weird RasterStack to have dims 
    # with no layers? but maybe it works?
    @testset "5.15 domain dimension" begin
        rast = RasterStack(examples["5.15"]; lazy=true)
        @test dims(rast) isa Tuple{<:Ti{<:Sampled{<:DateTime}},<:Dim{:pres},<:Y,<:X}
        @test layers(rast) == (;)
        @test refdims(rast) == ()
    end
    @testset "5.16 domain dimension" begin
        rast = RasterStack(examples["5.16"]; lazy=true)
        # Rotated pole dimensions are now typed as X and Y
        @test dims(rast) isa Tuple{<:Z,<:Y,<:X}
        @test layers(rast) == (;)
        @test refdims(rast) isa Tuple{<:Ti}
    end
    @testset "5.17 domain defines dimension and coords" begin
        # TODO : where to put area data? It could be useful
        # but currently has no formal location, so its just available
        # to the user, rather than to the `Rasters.cellarea` function
        # TODO and do we make polygons out of the vertices so this is a GeometryLookup?
        # Needs actual data in the file to test this
        rast = RasterStack(examples["5.17"]; lazy=true)
        span(dims(dims(rast, :cell), X))
    end
    @testset "5.17 domain with no layers or dimensions, but single coordinates, so refdims" begin
        rast = RasterStack(examples["5.18"]; lazy=true)
        @test dims(rast) == ()
        @test refdims(rast) == (DimensionalData.format(Ti(Sampled([0.5]; span=Regular(0.0)))),)
    end
    @testset "5.17 domain with timeseries geometries" begin
        rast = RasterStack(examples["5.19"]; lazy=true)
        @test lookup(rast) isa Tuple{<:GeometryLookup,Sampled}
        @test dims(rast) isa Tuple{<:Geometry,Ti}
        @test length(lookup(rast, Geometry)) == 2
        @test lookup(rast, Geometry)[1] == GeoInterface.LineString{false, false}([(30.0, 10.0), (10.0, 30.0), (40.0, 40.0)])
    end
end

@testset "5.20 indexed ragged array " begin
    # TODO ragged array lookup
    rast = RasterStack(examples["5.20"]; lazy=true)
    dims(rast, :station)
end

@testset "5.21 mesh" begin
    rast = RasterStack(examples["5.21"]; lazy=true)
end

@testset "taxon name/id" begin
    rast = Raster(examples["6.1.2"]; lazy=true)
    @test name(rast) == :abundance
    # Multiple char categorical coords are now merged into a MergedLookup
    taxon_lookup = lookup(rast, :taxon)
    @test taxon_lookup isa MergedLookup
    @test taxon_lookup[1] == ("urn:lsid:marinespecies.org:taxname:104464", "Calanus finmarchicus")
    # Inner dimensions accessible via .dims
    inner_dims = taxon_lookup.dims
    @test name(inner_dims[1]) == :taxon_lsid
    @test name(inner_dims[2]) == :taxon_name
    @test lookup(inner_dims[1]) isa Categorical{String,Vector{String},Unordered}
    @test lookup(inner_dims[1])[1] == "urn:lsid:marinespecies.org:taxname:104464"
    @test metadata(lookup(inner_dims[1]))["standard_name"] == "biological_taxon_lsid"
    @test lookup(inner_dims[2])[1] == "Calanus finmarchicus"
    # Selection by full tuple
    @test size(rast[taxon=At(("urn:lsid:marinespecies.org:taxname:104464", "Calanus finmarchicus"))]) == (100,)
    # Selection by inner dimension selectors
    @test size(rast[taxon=(Dim{:taxon_lsid}(At("urn:lsid:marinespecies.org:taxname:104464")), Dim{:taxon_name}(At("Calanus finmarchicus")))]) == (100,)
    # Selection by direct inner dimension name — same DiskArrays LogicalIndex bug as 5.3
    @test_broken (size(rast[Dim{:taxon_lsid}(At("urn:lsid:marinespecies.org:taxname:104464"))]) == (1, 100))
    @test_broken (size(rast[Dim{:taxon_name}(At("Calanus finmarchicus"))]) == (1, 100))
    # Selection by Where on tuple fields
    @test size(rast[taxon=Where(t -> occursin("finmarchicus", t[2]))]) == (1, 100)
end

@testset "6.1 region" begin
    rast = RasterStack(examples["6.1"]; lazy=true)
    @test name(rast) == (:n_heat_transport,)
    @test lookup(rast, :geo_region) isa Categorical{String,Vector{String},Unordered}
    @test lookup(rast, :geo_region)[1] == "atlantic_ocean"
    @test metadata(lookup(rast, :geo_region))["standard_name"] == "region"
    @test lookup(rast, Y) == 10.0:10:50
    @test lookup(rast, Ti)[1] == DateTime(1990)
end

@testset "6.2 secondary lookup coordinates" begin
    rast = Raster(examples["6.2"]; lazy=true)
    @test name(rast) == :xwind
    # Dimension variable + auxiliary coordinate are now merged into a MergedLookup
    sigma_lookup = lookup(rast, :sigma)
    @test sigma_lookup isa MergedLookup
    @test sigma_lookup[1] == (-1.0f0, Int32(10))  # (sigma value, model_level value)
    # Inner dimensions accessible via .dims
    inner_dims = sigma_lookup.dims
    @test name(inner_dims[1]) == :sigma
    @test name(inner_dims[2]) == :model_level
    @test collect(lookup(inner_dims[1])) == [-1.0f0, -2.0f0, -3.0f0, -4.0f0, -5.0f0]
    @test collect(lookup(inner_dims[2])) == Int32[10, 20, 30, 40, 50]
    # Selection by full tuple (model_level is Int32)
    @test size(rast[sigma=At((-1.0f0, Int32(10)))]) == (30,)
    # Selection by inner dimension selectors
    @test size(rast[sigma=(Dim{:sigma}(At(-1.0f0)), Dim{:model_level}(At(Int32(10))))]) == (30,)
    # Selection by direct inner dimension name (model_level) - dims are (Y, sigma)
    # Same DiskArrays LogicalIndex bug as 5.3 / taxon name/id
    @test_broken (size(rast[Dim{:model_level}(At(Int32(10)))]) == (30, 1))
    @test_broken (size(rast[Dim{:model_level}(At(Int32(30)))]) == (30, 1))
    # Selection by Where on model_level values
    @test size(rast[sigma=Where(t -> t[2] >= 30)]) == (30, 3)
end

@testset "7.2 non-aligned horizontal grid with polygon bounds" begin
    rast = RasterStack(examples["7.2"])
    @test haskey(rast, :dat)

    dat = rast[:dat]
    # 2D grid structure is preserved
    @test size(dat) == (2, 4)
    @test hasdim(dat, Dim{:imax})
    @test hasdim(dat, Dim{:jmax})

    # Check the ProjectedArrayLookup has geom_lookup field
    l1 = lookup(dims(dat, 1))
    @test l1 isa Rasters.ProjectedArrayLookup
    @test !isnothing(l1.geom_lookup)
    @test l1.geom_lookup isa GeometryLookup
    @test length(l1.geom_lookup) == 8  # 2x4 grid = 8 cells

    # Test Contains selector works (requires both dimensions)
    # Cell (1,1) has lat bounds 5-15, lon bounds 95-105, value 1.0
    point1 = GeoInterface.Point((100.0, 10.0))  # Inside cell (1,1)
    result1 = dat[Dim{:imax}(Contains(point1)), Dim{:jmax}(Contains(point1))]
    @test result1 == Float32(1.0)

    # Cell (2,1): lon bounds 105-115, value 2.0
    point2 = GeoInterface.Point((110.0, 10.0))  # Inside cell (2,1)
    result2 = dat[Dim{:imax}(Contains(point2)), Dim{:jmax}(Contains(point2))]
    @test result2 == Float32(2.0)

    # Normal integer indexing still works
    @test dat[1, 1] == Float32(1.0)
    @test dat[2, 1] == Float32(2.0)
end

@testset "7.3 formula terms" begin
    # These are not implemented, but they load ok and warn
    @test_warn "formula_terms" RasterStack(examples["7.3"]; lazy=true)
end

@testset "7.3 cell areas" begin
    # TODO load bounds as polygons
    # Not sure if we should attach the area to the dimension or its fine as a variable?
    RasterStack(examples["7.4"]; lazy=true)
end

@testset "7.5 methods applied to a timeseries" begin
    rast = RasterStack(examples["7.5"]; lazy=true)
    # TODO: unfortunately mixed points and intervals 
    # on the same axis is hard to represent in DD without a DimTree
    @test metadata(rast.ppn)["cell_methods"] == "time: sum"
    @test metadata(rast.pressure)["cell_methods"] == "time: point"
    @test metadata(rast.maxtemp)["cell_methods"] == "time: maximum"
end

@testset "7.6 spacing of data" begin
    rast = RasterStack(examples["7.6"]; lazy=true)
    @test lookup(rast, Ti) == [DateTime(1990, 1, 1, 12)]
    @test Rasters.bounds(rast, Ti) == (DateTime(1990, 1, 1, 0), DateTime(1990, 1, 2, 0))
    @test metadata(rast.TS_var)["cell_methods"] == "time: variance (interval: 1 hr comment: sampled instantaneously)"
end

@testset "7.7 labelled categorical axis" begin
    rast = RasterStack(examples["7.7"]; lazy=true)
    @test lookup(rast, :land_sea) == ["land", "sea"]
end

@testset "climatology bounds" begin
    rast = RasterStack(examples["7.9"]; lazy=true)
    @test Rasters.Lookups.cycle(lookup(rast, Ti)) == Year(1)
    @test Rasters.bounds(rast, Ti) == (DateTime("1960-03-01T00:00:00"), DateTime("1991-03-01T00:00:00"))
    @test span(rast, Ti) == Explicit([
        DateTime("1960-03-01T00:00:00") DateTime("1960-06-01T00:00:00") DateTime("1960-09-01T00:00:00") DateTime("1960-12-01T00:00:00")
        DateTime("1960-06-01T00:00:00") DateTime("1960-09-01T00:00:00") DateTime("1960-12-01T00:00:00") DateTime("1961-03-01T00:00:00")
    ])
    # Near selects the closest interval - March is closest to the first interval (March-June)
    @test rast[Ti=Near(DateTime(1960, 3, 15))] == rast[Ti=1]
    @test rast[Ti=At(DateTime(1960, 4, 16))] ==
          rast[Ti=At(DateTime(1970, 4, 16))] ==
          rast[Ti=At(DateTime(1990, 4, 16))] == rast[Ti=1]
    # Years outside climatology bounds (1960-1991) should error, not wrap
    @test try rast[Ti=At(DateTime(2000, 4, 16))]; false catch e; e isa Lookups.SelectorError end
    @test try rast[Ti=At(DateTime(1950, 4, 16))]; false catch e; e isa Lookups.SelectorError end
    # Contains is just broken for Cyclic, this should work
    @test rast[Ti=Contains(DateTime(1970, 6, 1))] == rast[Ti=2]
    @test rast[Ti=Contains(DateTime(1970, 5, 31))] == rast[Ti=1]

    rast = RasterStack(examples["7.10"]; lazy=false)
    @test Rasters.Lookups.cycle(lookup(rast, Ti)) == Year(1)
    @test Rasters.bounds(rast, Ti) == (DateTime("1961-01-01T00:00:00"), DateTime("1990-02-01T00:00:00"))
    @test span(rast, Ti) == Explicit([
        DateTime("1961-01-01T00:00:00") DateTime("1971-01-01T00:00:00") DateTime("1981-01-01T00:00:00")
        DateTime("1961-02-01T00:00:00") DateTime("1971-02-01T00:00:00") DateTime("1981-02-01T00:00:00")
    ])
    # Mix up the array values
    @test Rasters.dims2indices(rast, (Ti(At(DateTime(1961, 1, 15))),)) == (:, :, 1)
    @test_broken Rasters.dims2indices(rast, (Ti(At(DateTime(1971, 1, 15))),)) == (:, :, 2)
    @test_broken Rasters.dims2indices(rast, (Ti(At(DateTime(1981, 1, 15))),)) == (:, :, 3)
    rast = RasterStack(examples["7.11"]; lazy=true)
    @test Rasters.Lookups.cycle(lookup(rast, Ti)) == Day(1)
    @test Rasters.bounds(rast, Ti) == (DateTime("1997-04-01T00:00:00"), DateTime("1997-05-02T00:00:00"))
    rast = RasterStack(examples["7.12"]; lazy=true)
    @test Rasters.Lookups.cycle(lookup(rast, Ti)) == Year(1)
    @test Rasters.bounds(rast, Ti) == (DateTime("2007-12-01T06:00:00"), DateTime("2008-3-01T06:00:00"))
    rast = RasterStack(examples["7.13"]; lazy=true)
    @test Rasters.Lookups.cycle(lookup(rast, Ti)) == (Year(1), Day(1) => Month(1))
    @test Rasters.bounds(rast, Ti) == (DateTime("1961-04-01T00:00:00"), DateTime("1990-05-01T00:00:00"))
    rast = RasterStack(examples["7.14"]; lazy=true)
    @test Rasters.bounds(rast, Ti) == (DateTime("2000-06-01T00:00:00"), DateTime("2000-09-01T00:00:00"))
    @test Rasters.Lookups.cycle(lookup(rast, Ti)) == Year(1)
end

@testset "Geometry lookups" begin
    @testset "7.15 linestring" begin
        rast = RasterStack(examples["7.15"])
        @test size(rast) == (4, 2)
        @test dims(rast, :Geometry)[1] isa GeoInterface.Wrappers.LineString
        @test GeoInterface.getpoint(dims(rast, :Geometry)[1]) == [(30.0, 10.0), (10.0, 30.0), (40.0, 40.0)]
    end
    @testset "7.16 polygons with holes" begin
        rast = RasterStack(examples["7.16"])
        multipoly = dims(rast, :Geometry)[1]
        @test multipoly isa GeoInterface.MultiPolygon
        poly = GeoInterface.getgeom(dims(rast, :Geometry)[1], 1)
        @test poly isa GeoInterface.Polygon
        GeoInterface.nring(poly) == 2
        GeoInterface.nhole(poly) == 1
        @test GeoInterface.getgeom(dims(rast, :Geometry)[1], 1) isa GeoInterface.Polygon
        @test collect(GeoInterface.getpoint(dims(rast, :Geometry)[1])) ==
            [(20.0, 0.0), (10.0, 15.0), (0.0, 0.0), (20.0, 0.0),
             (5.0, 5.0), (10.0, 10.0), (15.0, 5.0), (5.0, 5.0),
             (20.0, 20.0), (10.0, 35.0), (0.0, 20.0), (20.0, 20.0)]
    end
end

@testset "CRS round-trip" begin
    # Test writing CRS to NetCDF and reading it back
    data = rand(Float32, 10, 10)
    rast = Raster(data, (X(1.0:10.0), Y(1.0:10.0)); crs=EPSG(32615), name=:temperature)

    mktempdir() do dir
        path = joinpath(dir, "crs_roundtrip_test.nc")
        write(path, rast)
        rast2 = Raster(path)
        @test !isnothing(crs(rast2))
        # CRS should be preserved (converted to WKT and back)
        @test crs(rast2) isa WellKnownText2
    end
end

# ---------------------------------------------------------------------------
# DimensionalData syntax access:
#
# Goal: validate the *user-facing* DD experience for every loaded CF example.
# A CF NetCDF should drop into Rasters/DimensionalData as a coherent typed
# stack with named dimensions, working selectors (At, Near, Contains, Where,
# ..), and correct dim reduction.
#
# For examples where the CF spec stored real values (e.g. 5.1, 5.6, 5.10, 7.2,
# 7.15, 7.16) we check the *values* round-tripped through the selector are
# correct. For the structural CF examples that ship filled with the CF default
# fill value (9.96921f36) we only check the access pattern works and returns
# the right shape/type/dim reduction.
# ---------------------------------------------------------------------------

@testset "DimensionalData access" begin

    @testset "2.1 string variables" begin
        rast = RasterStack(examples["2.1"]; lazy=true)
        @test ndims(rast.str_variable) == 1
        @test size(rast.str_variable) == (30,)
        # NoLookup so only integer indexing
        @test rast.str_variable[1] isa Union{Missing,AbstractString}
        @test size(rast.str_variable[1:10]) == (10,)
        # Stack key access
        @test rast[:str_variable] === rast.str_variable
    end

    @testset "3.1 zero-dim scalar layers" begin
        rast = RasterStack(examples["3.1"])
        @test ndims(rast.Tonscale) == 0
        @test ndims(rast.Tdifference) == 0
        @test rast.Tonscale[] isa Union{Missing,Float32}
    end

    @testset "5.1 independent coords - full selector matrix" begin
        rast = RasterStack(examples["5.1"])
        x = rast.xwind
        @test size(x) == (36, 18, 15, 4)

        # Integer indexing
        @test ndims(x[1, 1, 1, 1]) == 0
        @test size(x[1:10, :, :, :]) == (10, 18, 15, 4)

        # Named-dim integer
        @test size(x[X=1]) == (18, 15, 4)
        @test !hasdim(x[X=1], X)

        # At selector reduces dim
        sub = x[X=At(1.0f0)]
        @test size(sub) == (18, 15, 4)
        @test !hasdim(sub, X)
        @test hasdim(sub, Y) && hasdim(sub, Dim{:pres}) && hasdim(sub, Ti)

        # Near picks closest
        @test x[X=Near(5.4f0)] == x[X=5]
        @test x[Y=Near(10.0f0)] == x[Y=10]

        # Interval selector
        @test size(x[X=5.0f0..10.0f0]) == (6, 18, 15, 4)
        @test size(x[Y=1.0f0..5.0f0]) == (36, 5, 15, 4)

        # Ti selectors
        @test x[Ti=At(DateTime(1990,1,2))] == x[Ti=2]
        @test size(x[Ti=DateTime(1990,1,1)..DateTime(1990,1,3)]) == (36, 18, 15, 3)

        # Multiple combined selectors
        sub = x[X=At(1.0f0), Y=At(1.0f0), Ti=At(DateTime(1990,1,1))]
        @test size(sub) == (15,)
        @test hasdim(sub, Dim{:pres})

        # Pres dim by name (Dim{:pres})
        @test size(x[pres=At(1.0f0)]) == (36, 18, 4)
        @test size(x[pres=1..5]) == (36, 18, 5, 4)
    end

    @testset "5.2 2D coordinate variables (ProjectedArrayLookup)" begin
        rast = RasterStack(examples["5.2"]; lazy=true)
        T = rast.T
        # Data is fill - test structural access only
        @test size(T) == (128, 64, 18)
        @test size(T[X=1, Y=1]) == (18,)
        @test size(T[X=1]) == (64, 18)
        @test size(T[Y=1]) == (128, 18)
        @test size(T[Z=1]) == (128, 64)
        @test !hasdim(T[Z=1], Z)
        @test size(T[Z=1:5]) == (128, 64, 5)
    end

    @testset "5.3 reduced grid (MergedLookup)" begin
        rast = RasterStack(examples["5.3"]; lazy=true)
        @test rast.PS[rgrid=1] === 0.1f0
        @test rast.PS[rgrid=4] === 0.4f0
        # Selector on full tuple
        @test rast.PS[rgrid=At((1.0, 10.0))] === 0.1f0
        @test rast.PS[rgrid=At((4.0, 40.0))] === 0.4f0
        # Where on inner-tuple fields
        @test size(rast.PS[rgrid=Where(t -> t[1] >= 3.0)]) == (2,)
        @test collect(rast.PS[rgrid=Where(t -> t[2] == 30.0)]) == Float32[0.3]
    end

    @testset "5.6 rotated pole grid (ProjectedArrayLookup)" begin
        rast = RasterStack(examples["5.6"]; lazy=true)
        t = rast.temp
        # Integer access
        @test t[X=1, Y=1, Z=1] === 1.0f0
        @test t[X=2, Y=1, Z=1] === 2.0f0
        @test t[X=1, Y=2, Z=1] === 3.0f0
        @test t[X=2, Y=3, Z=2] === 12.0f0
        # Z selectors return correct shape and value
        @test size(t[Z=At(100.0f0)]) == (2, 3)
        @test t[Z=At(100.0f0)] == t[Z=1]
        @test t[Z=Near(140.0f0)] == t[Z=1]
        # Geographic Near via PAL matrix
        @test t[X=Near(-0.9), Y=Near(3.0), Z=1] === 3.0f0
    end

    @testset "5.8-9 simple lat/lon (fill data)" begin
        for id in ("5.8", "5.9")
            rast = RasterStack(examples[id]; lazy=true)
            t = rast.temp
            @test size(t) == (36, 18)
            @test size(t[X=1, Y=1]) == ()
            @test size(t[X=1:10]) == (10, 18)
        end
    end

    @testset "5.10 British national grid - integer + Z selectors" begin
        rast = RasterStack(examples["5.10"]; lazy=true)
        t = rast.temp
        @test t[X=1, Y=1, Z=1] === 1.0f0
        @test t[X=2, Y=3, Z=2] === 12.0f0
        @test size(t[Z=At(100.0f0)]) == (2, 3)
        @test size(t[Z=1]) == (2, 3)
        # mappedcrs accessible per-dim
        @test !isnothing(lookup(t, X).mappedcrs)
    end

    @testset "5.11 WGS84 - real X/Y selectors" begin
        rast = RasterStack(examples["5.11"]; lazy=true)
        t = rast.temp
        @test size(t) == (36, 18)
        # X is 0..350 step 10, Y is -85..85 step 10
        @test size(t[X=At(0.0)]) == (18,)
        @test size(t[Y=At(-85.0)]) == (36,)
        @test size(t[X=Near(5.0)]) == (18,)
        # Interval: X is 0,10,...,350 → 0..100 = 11 values
        @test size(t[X=0.0..100.0]) == (11, 18)
        # Y is -85,-75,...,85 → -50..50 grabs -45..45 = 10 values (50 isn't on grid)
        @test size(t[Y=-50.0..50.0]) == (36, 10)
        # Combined
        @test size(t[X=0.0..100.0, Y=-50.0..50.0]) == (11, 10)
    end

    @testset "5.14 refdims accessible" begin
        rast = RasterStack(examples["5.14"]; lazy=true)
        @test length(refdims(rast)) == 2
        h = rast.height
        # Ti selector reduces dim and adds to refdims
        sub = h[Ti=At(DateTime(1999, 1, 1, 12))]
        @test !hasdim(sub, Ti)
        @test hasdim(sub, X) && hasdim(sub, Y)
        @test size(sub) == (360, 180)
        # Ti Near
        @test h[Ti=Near(DateTime(1999, 1, 1, 12, 30))] == h[Ti=2]
        # Ti interval
        @test size(h[Ti=DateTime(1999, 1, 1, 6)..DateTime(1999, 1, 1, 18)]) == (360, 180, 3)
    end

    @testset "5.15 domain stack - dims without layers" begin
        rast = RasterStack(examples["5.15"]; lazy=true)
        @test layers(rast) == (;)
        @test length(dims(rast)) == 4
        @test hasdim(rast, X) && hasdim(rast, Y) && hasdim(rast, Ti) && hasdim(rast, Dim{:pres})
        # Ti values are real
        @test lookup(rast, Ti)[1] == DateTime(1990, 1, 2)
    end

    @testset "5.16 domain stack - rotated pole, refdims" begin
        rast = RasterStack(examples["5.16"]; lazy=true)
        @test layers(rast) == (;)
        @test hasdim(rast, X) && hasdim(rast, Y) && hasdim(rast, Z)
        @test length(refdims(rast)) == 1
    end

    @testset "5.17 cell areas with MergedLookup cell + Ti" begin
        rast = RasterStack(examples["5.17"]; lazy=true)
        ca = rast.cell_area
        # cell_area is 1D (just :cell), not over Ti
        @test size(ca) == (2562,)
        @test size(ca[cell=1:10]) == (10,)
        @test ndims(ca[cell=1]) == 0
        # MergedLookup tuple data
        @test lookup(rast, :cell) isa MergedLookup
    end

    @testset "5.18 refdims only, no dims" begin
        rast = RasterStack(examples["5.18"]; lazy=true)
        @test dims(rast) == ()
        @test length(refdims(rast)) == 1
        @test name(refdims(rast)[1]) === :Ti
    end

    @testset "5.19 timeseries geometry" begin
        rast = RasterStack(examples["5.19"]; lazy=true)
        @test layers(rast) == (;)
        @test hasdim(rast, :Geometry) && hasdim(rast, Ti)
        @test length(lookup(rast, :Geometry)) == 2
        @test lookup(rast, :Geometry)[1] == GeoInterface.LineString{false,false}(
            [(30.0, 10.0), (10.0, 30.0), (40.0, 40.0)]
        )
    end

    @testset "5.20 indexed ragged - station NoLookup" begin
        rast = RasterStack(examples["5.20"]; lazy=true)
        @test size(rast.station_info) == (23,)
        @test size(rast.station_info[1:5]) == (5,)
        @test rast.station_info[station=1] isa Union{Missing,Int32}
    end

    @testset "5.21 UGRID mesh" begin
        rast = RasterStack(examples["5.21"]; lazy=true)
        # Each layer has different dim shape
        @test size(rast.mesh_face_nodes) == (4, 2)
        @test size(rast.mesh_edge_nodes) == (2, 6)
        # Ti selectors on layers that have Ti
        @test size(rast.volume_at_faces) == (2, 12)
        @test size(rast.volume_at_faces[Ti=1]) == (2,)
        @test size(rast.volume_at_faces[Ti=At(DateTime(2004, 6, 5))]) == (2,)
        @test size(rast.volume_at_faces[Ti=DateTime(2004, 6, 5)..DateTime(2004, 6, 9)]) == (2, 5)
    end

    @testset "6.1 region categorical Y selectors" begin
        rast = RasterStack(examples["6.1"]; lazy=true)
        nht = rast.n_heat_transport
        @test size(nht) == (1, 5, 4)
        # Categorical selector
        @test size(nht[geo_region=At("atlantic_ocean")]) == (5, 4)
        @test !hasdim(nht[geo_region=At("atlantic_ocean")], :geo_region)
        # Y selector (Mapped float)
        @test size(nht[Y=At(20.0f0)]) == (1, 4)
        @test size(nht[Y=Near(25.0f0)]) == (1, 4)
        # Y interval
        @test size(nht[Y=20.0f0..40.0f0]) == (1, 3, 4)
        # Ti selector
        @test size(nht[Ti=At(DateTime(1990, 1, 1))]) == (1, 5)
    end

    @testset "6.1.2 taxon - full tuple + inner-dim selectors" begin
        rast = Raster(examples["6.1.2"]; lazy=true)
        # Selection by full tuple
        sub = rast[taxon=At(("urn:lsid:marinespecies.org:taxname:104464", "Calanus finmarchicus"))]
        @test size(sub) == (100,)
        @test hasdim(sub, Ti)
        # Inner-dim tuple
        sub2 = rast[taxon=(
            Dim{:taxon_lsid}(At("urn:lsid:marinespecies.org:taxname:104464")),
            Dim{:taxon_name}(At("Calanus finmarchicus")),
        )]
        @test size(sub2) == (100,)
        # Where on inner fields
        @test size(rast[taxon=Where(t -> occursin("finmarchicus", t[2]))]) == (1, 100)
        # Ti selector
        @test size(rast[Ti=At(DateTime(2019, 1, 5))]) == (2,)
        @test size(rast[Ti=Near(DateTime(2019, 2, 1))]) == (2,)
        # Ti interval
        @test size(rast[Ti=DateTime(2019, 1, 2)..DateTime(2019, 1, 10)]) == (2, 9)
    end

    @testset "6.2 model_level MergedLookup" begin
        rast = Raster(examples["6.2"]; lazy=true)
        @test size(rast) == (30, 5)
        # Full tuple
        sub = rast[sigma=At((-1.0f0, Int32(10)))]
        @test size(sub) == (30,)
        # Inner-dim tuple
        sub2 = rast[sigma=(
            Dim{:sigma}(At(-1.0f0)),
            Dim{:model_level}(At(Int32(10))),
        )]
        @test size(sub2) == (30,)
        # Where on second tuple field
        @test size(rast[sigma=Where(t -> t[2] >= 30)]) == (30, 3)
        # Y is NoLookup so integer only
        @test size(rast[Y=1]) == (5,)
        @test size(rast[Y=1:5]) == (5, 5)
    end

    @testset "7.2 non-aligned grid - geographic Contains via bounds" begin
        rast = RasterStack(examples["7.2"])
        d = rast.dat
        @test size(d) == (2, 4)
        # Integer access matches values
        @test d[1, 1] === 1.0f0
        @test d[2, 1] === 2.0f0
        @test d[1, 4] === 7.0f0
        @test d[2, 4] === 8.0f0
        # Geographic Contains via geom_lookup
        point1 = GeoInterface.Point((100.0, 10.0))
        @test d[Dim{:imax}(Contains(point1)), Dim{:jmax}(Contains(point1))] === 1.0f0
        point2 = GeoInterface.Point((110.0, 10.0))
        @test d[Dim{:imax}(Contains(point2)), Dim{:jmax}(Contains(point2))] === 2.0f0
        # Subrange slicing
        @test size(d[imax=1:2, jmax=1:2]) == (2, 2)
    end

    @testset "7.3 formula terms" begin
        rast = RasterStack(examples["7.3"]; lazy=true)
        # eta is MergedLookup; X/Y NoLookup
        @test size(rast.temp) == (10, 10, 5)
        @test size(rast.temp[X=1, Y=1]) == (5,)
        @test size(rast.temp[eta=1]) == (10, 10)
        @test size(rast.PS) == (10, 10)
    end

    @testset "7.4 cell areas with Ti" begin
        rast = RasterStack(examples["7.4"]; lazy=true)
        @test size(rast.PS) == (2562, 12)
        @test size(rast.PS[Ti=At(DateTime(1979, 1, 2))]) == (2562,)
        @test size(rast.PS[Ti=Near(DateTime(1979, 1, 7, 12))]) == (2562,)
        @test size(rast.PS[Ti=DateTime(1979, 1, 2)..DateTime(1979, 1, 5)]) == (2562, 4)
        # cell_area has no Ti
        @test size(rast.cell_area) == (2562,)
        @test size(rast.cell_area[cell=1:100]) == (100,)
    end

    @testset "7.5 timeseries methods" begin
        rast = RasterStack(examples["7.5"]; lazy=true)
        @test size(rast.pressure) == (10, 5)
        @test size(rast.pressure[Ti=At(DateTime(1998, 4, 19, 6))]) == (10,)
        @test size(rast.pressure[station=1]) == (5,)
        # Ti is 5 stamps step 12h. Interval bounds are half-open at the upper
        # edge (Intervals sampling) so 06:00..18:00 +1day grabs 3 of 5.
        @test size(rast.pressure[Ti=DateTime(1998, 4, 19, 6)..DateTime(1998, 4, 20, 18)]) == (10, 3)
        # Last stamp via At
        @test size(rast.pressure[Ti=At(DateTime(1998, 4, 21, 6))]) == (10,)
    end

    @testset "7.6 single Ti" begin
        rast = RasterStack(examples["7.6"]; lazy=true)
        v = rast.TS_var
        @test size(v) == (180, 90, 1)
        @test size(v[Ti=1]) == (180, 90)
        @test size(v[Ti=At(DateTime(1990, 1, 1, 12))]) == (180, 90)
        # Contains based on bounds (00:00..24:00)
        @test size(v[Ti=Contains(DateTime(1990, 1, 1, 6))]) == (180, 90)
    end

    @testset "7.7 categorical land_sea" begin
        rast = RasterStack(examples["7.7"]; lazy=true)
        flux = rast.surface_upward_sensible_heat_flux
        @test size(flux) == (96, 73, 2)
        @test size(flux[land_sea=At("land")]) == (96, 73)
        @test size(flux[land_sea=At("sea")]) == (96, 73)
        @test !hasdim(flux[land_sea=At("land")], :land_sea)
        # surface_temperature does not have land_sea
        @test size(rast.surface_temperature) == (96, 73)
    end

    @testset "7.9-7.14 climatology Cyclic Ti" begin
        # 7.9 seasonal (already deeply tested) - add cyclic year wrap probe
        rast = RasterStack(examples["7.9"]; lazy=true)
        # In-bounds: any year wraps
        @test rast[Ti=At(DateTime(1960, 4, 16))] == rast[Ti=At(DateTime(1985, 4, 16))]
        # Contains works
        @test rast[Ti=Contains(DateTime(1970, 6, 1))] == rast[Ti=2]

        # 7.11 hourly cycle - X/Y are NoLookup so integer only
        rast = RasterStack(examples["7.11"]; lazy=true)
        t = rast.temperature
        @test size(t) == (10, 10, 24)
        @test size(t[Ti=At(DateTime(1997, 4, 1, 12, 30))]) == (10, 10)
        # Cycle is Day(1): same hour-of-day on a different in-bounds date
        # should return the same slice
        @test t[Ti=At(DateTime(1997, 4, 1, 12, 30))] ==
              t[Ti=At(DateTime(1997, 4, 15, 12, 30))]
        # A climatology is bounded to its reference period. Querying outside
        # those bounds (e.g. 2025) correctly errors - extrapolating a 1997
        # climatology to 2025 would silently assume climate stability.
        @test_throws Lookups.SelectorError t[Ti=At(DateTime(2025, 4, 1, 12, 30))]

        # 7.10 decadal - integer indexing maps to right slice
        rast = RasterStack(examples["7.10"]; lazy=true)
        @test Rasters.dims2indices(rast, (Ti(At(DateTime(1961, 1, 15))),)) == (:, :, 1)
    end

    @testset "7.15 LineString geometries x Ti" begin
        rast = RasterStack(examples["7.15"])
        sd = rast.someData
        @test size(sd) == (4, 2)
        # Full selection (every dim specified) works
        @test sd[Ti=1, Geometry=1] == 1.0
        @test sd[Ti=1, Geometry=2] == 1.0
        @test sd[Ti=At(DateTime(2000, 1, 4)), Geometry=2] == 3.0
        # Partial selection (only Ti or only Geometry) works
        @test size(sd[Geometry=1]) == (4,)
        @test size(sd[Ti=1]) == (2,)
        @test size(sd[Ti=DateTime(2000, 1, 2)..DateTime(2000, 1, 3)]) == (2, 2)
        # Geometry lookup is GeoInterface objects
        @test dims(rast, :Geometry)[1] isa GeoInterface.Wrappers.LineString
    end

    @testset "7.16 MultiPolygon geometries x Ti" begin
        rast = RasterStack(examples["7.16"])
        sd = rast.someData
        @test size(sd) == (4, 2)
        @test sd[Ti=1, Geometry=1] == 1.0
        @test sd[Ti=At(DateTime(2000, 1, 4)), Geometry=2] == 3.0
        @test dims(rast, :Geometry)[1] isa GeoInterface.MultiPolygon
    end

end

# ---------------------------------------------------------------------------
# Round-trip every CF example: load -> write -> load -> compare.
#
# Compares structure (dim names, lookup type families, layer keys) and data.
# CRS encoding is allowed to change (EPSG -> WKT), so only its presence is
# checked, not its exact type. Per-example assertions are tagged @test or
# @test_broken depending on whether the current write path supports the
# example - this surfaces gaps in CI instead of hiding them.
# ---------------------------------------------------------------------------

# Compare two lookups loosely: lookup-type family and values match.
# Returns true if compatible; raises @test failure with useful message otherwise.
function _lookups_match(la, lb)
    typeof(la).name === typeof(lb).name || return false
    # Values can differ in eltype after round-trip (e.g. Float32 vs Float64)
    # but should compare ==
    try
        collect(la) == collect(lb)
    catch
        false
    end
end

function _dims_match(da, db)
    length(da) == length(db) || return false
    for (a, b) in zip(da, db)
        typeof(a).name === typeof(b).name || return false
        name(a) === name(b) || return false
        _lookups_match(lookup(a), lookup(b)) || return false
    end
    return true
end

function _layers_match(sa, sb)
    keys(sa) == keys(sb) || return false
    for k in keys(sa)
        a = read(sa[k])
        b = read(sb[k])
        size(a) == size(b) || return false
        # NaN-aware equality
        all(((x, y),) -> isequal(x, y), zip(parent(a), parent(b))) || return false
    end
    return true
end

function _roundtrip(path, tmpdir)
    out = joinpath(tmpdir, "rt_" * basename(path))
    orig = RasterStack(path; lazy=false)
    write(out, orig; force=true)
    reloaded = RasterStack(out; lazy=false)
    return orig, reloaded, out
end

# ---- dangling-reference check ----
# Walks every CF reference-style attribute in a written file and confirms the
# variable it points at actually exists. Catches silent gaps where the writer
# preserves a `bounds=lon_vertices` attribute but never writes `lon_vertices`.
function _parse_cf_refs(key::AbstractString, val)
    val isa AbstractString || return String[]
    if key in ("bounds", "climatology", "geometry", "nodes",
               "node_count", "node_coordinates", "part_node_count",
               "interior_ring")
        # Some of these are space-separated, some single. Split & take all.
        return String.(split(val))
    elseif key in ("coordinates", "ancillary_variables")
        return String.(split(val))
    elseif key == "grid_mapping"
        # "name" or "name1: c1 c2 name2: c3 c4" — pull the names before ':'
        if occursin(":", val)
            return [String(rstrip(t, ':')) for t in split(val) if endswith(t, ":")]
        else
            return [String(val)]
        end
    elseif key in ("cell_measures", "formula_terms")
        # "k1: v1 k2: v2" — take the v's (even-indexed tokens)
        toks = split(val)
        return [String(t) for (i, t) in enumerate(toks) if iseven(i) && !endswith(t, ":")]
    end
    return String[]
end

function _dangling_refs(path::AbstractString)
    out = Tuple{String,String,String}[]
    NCDatasets.NCDataset(path) do ds
        for varname in keys(ds)
            attrs = Dict(ds[varname].attrib)
            for (k, v) in attrs
                for ref in _parse_cf_refs(string(k), v)
                    haskey(ds, ref) || push!(out, (varname, string(k), ref))
                end
            end
        end
    end
    return out
end

# ---- metadata round-trip helpers ----
# Compare two metadata sources as dicts. `ignore` is a set of keys that are
# legitimately synthesized on write or stripped on read - their presence /
# value is not enforced. Any other key must match exactly in both directions.
_md_dict(::Lookups.NoMetadata) = Dict{String,Any}()
_md_dict(m::Lookups.Metadata) = Dict{String,Any}(string(k) => v for (k, v) in m)
_md_dict(d::AbstractDict) = Dict{String,Any}(string(k) => v for (k, v) in d)

function _metadata_match(orig, reloaded; ignore=())
    mo, mr = _md_dict(orig), _md_dict(reloaded)
    for k in keys(mo)
        k in ignore && continue
        haskey(mr, k) || return false
        isequal(mo[k], mr[k]) || return false
    end
    for k in keys(mr)
        k in ignore && continue
        haskey(mo, k) || return false
    end
    return true
end

# Layer attribs that are synthesized on write or stripped on read.
# `_FillValue` is written from the missingval pair; `coordinates` is built
# from aux-coord names; `grid_mapping` from the CRS write; scale/offset are
# handled via the Mod and won't survive verbatim.
const _LAYER_IGNORE = Set(("_FillValue", "missing_value", "coordinates",
                          "grid_mapping", "scale_factor", "add_offset"))
# Dim attribs that are synthesized on write. `axis` is added by
# `_cdm_set_axis_attrib!`. `standard_name` is added as a default for Z and
# Ti dims when not already present (now uses `get!`, so values from the
# source are preserved) - ignore to avoid false positives where a source
# file omitted it.
const _DIM_IGNORE = Set(("_FillValue", "missing_value", "axis", "standard_name"))

function _layer_metadata_match(orig, reloaded)
    keys(orig) == keys(reloaded) || return false
    for k in keys(orig)
        _metadata_match(metadata(orig[k]), metadata(reloaded[k]);
                        ignore=_LAYER_IGNORE) || return false
    end
    return true
end

function _dim_metadata_match(orig, reloaded)
    for d in dims(orig)
        d_re = DimensionalData.dims(reloaded, DimensionalData.basetypeof(d))
        isnothing(d_re) && continue
        _metadata_match(metadata(d), metadata(d_re);
                        ignore=_DIM_IGNORE) || return false
        # Recurse into the inner dims that MergedLookup and ProjectedArrayLookup
        # hold - that is where the per-aux-coord attribs (bounds, long_name,
        # units etc.) live, and silently losing them is exactly the class of
        # bug this test is here to catch.
        _inner_dim_metadata_match(lookup(d), lookup(d_re)) || return false
    end
    return true
end

_inner_dim_metadata_match(la, lb) = true
function _inner_dim_metadata_match(la::MergedLookup, lb::MergedLookup)
    length(la.dims) == length(lb.dims) || return false
    for (a, b) in zip(la.dims, lb.dims)
        _metadata_match(metadata(a), metadata(b); ignore=_DIM_IGNORE) || return false
    end
    return true
end
function _inner_dim_metadata_match(la::Rasters.ProjectedArrayLookup,
                                   lb::Rasters.ProjectedArrayLookup)
    length(la.dims) == length(lb.dims) || return false
    for (a, b) in zip(la.dims, lb.dims)
        # The NearestNeighbors extension strips inner-dim metadata via
        # `format_unaligned`. When that happens we can't usefully compare.
        (a isa DimensionalData.Dimension{<:Lookups.Lookup} &&
         b isa DimensionalData.Dimension{<:Lookups.Lookup}) || continue
        _metadata_match(metadata(a), metadata(b); ignore=_DIM_IGNORE) || return false
    end
    return true
end

@testset "CF round-trip" begin
    # Per-example expectations: (write_ok, dims_ok, data_ok).
    # Flip a flag to true once the underlying bug is fixed; the testset will
    # immediately fail with "unexpected pass" so the @test_broken can be
    # removed.
    #
    # Known gaps that drive the broken flags:
    # - ProjectedArrayLookup write: 5.2 5.6 5.10 round-trip cleanly now,
    #   but 7.2's dims are plain `Dim` (not X/Y), so write skips aux coords
    #   and the load itself drops one of the two matrices for non-X/Y dims.
    # - Domain stacks (no layers) drop dims on write because write is
    #   layer-driven (5.15 5.16 5.19)
    # - MergedLookup write: 5.3 6.1.2 6.2 7.3 7.4 now round-trip cleanly.
    #   5.17 is a domain stack (no layers) - same bucket as 5.15/5.16/5.19.
    # (Climatology Cyclic 7.9-7.14 now round-trip via _def_lookup_var! on
    #  AbstractCyclic; the on-disk bounds matrix is reconstructed by adding
    #  the duration back to upper bounds.)
    # meta_ok: dim & layer metadata round-trips after ignoring synthesized
    # keys (_FillValue, coordinates, grid_mapping, axis, standard_name -
    # see _LAYER_IGNORE and _DIM_IGNORE above). Real losses like 7.4 dropping
    # `bounds` polygon vars are caught by leaving `bounds` outside the ignore
    # list.
    expectations = Dict(
        "2.1"   => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "3.1"   => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "5.1"   => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "5.2"   => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "5.3"   => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "5.6"   => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "5.8"   => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "5.9"   => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "5.10"  => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "5.11"  => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "5.14"  => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "5.15"  => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "5.16"  => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "5.17"  => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "5.18"  => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "5.19"  => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "5.20"  => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "5.21"  => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "6.1"   => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "6.1.2" => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "6.2"   => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "7.2"   => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "7.3"   => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "7.4"   => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "7.5"   => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "7.6"   => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "7.7"   => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "7.9"   => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "7.10"  => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "7.11"  => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "7.12"  => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "7.13"  => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
        "7.14"  => (write_ok=true,  dims_ok=true,  data_ok=true,  meta_ok=true, refs_ok=true),
    )
    mktempdir() do tmpdir
        for (id, exp) in pairs(expectations)
            @testset "$id" begin
                if !exp.write_ok
                    @test_broken (_roundtrip(examples[id], tmpdir); true)
                    continue
                end
                orig, reloaded, outpath = _roundtrip(examples[id], tmpdir)
                @test keys(reloaded) == keys(orig)
                if exp.dims_ok
                    @test _dims_match(dims(orig), dims(reloaded))
                else
                    @test_broken _dims_match(dims(orig), dims(reloaded))
                end
                # Only attempt a data comparison when shapes line up. When a
                # dim is dropped, layers reshape and pairwise comparison would
                # throw - which would mask the underlying dim bug.
                shapes_match = keys(reloaded) == keys(orig) &&
                    all(k -> size(orig[k]) == size(reloaded[k]), keys(orig))
                if exp.data_ok && shapes_match
                    @test _layers_match(orig, reloaded)
                else
                    @test_broken (shapes_match && _layers_match(orig, reloaded))
                end
                # CRS presence is preserved when both dims write correctly
                # (it travels with grid_mapping on the spatial dims/layers).
                # When the dim write path is broken, the CRS is lost too —
                # mark broken in that case so the gap stays visible.
                # Multi-CRS (Vector{<:Pair}, e.g. CF 5.10 with both an OSGB
                # projection and WGS84 lat/lon) isn't supported by the writer
                # yet, so flag it broken instead of asserting.
                if !isnothing(crs(orig))
                    if crs(orig) isa AbstractVector
                        @test_broken !isnothing(crs(reloaded))
                    elseif exp.dims_ok
                        @test !isnothing(crs(reloaded))
                    else
                        @test_broken !isnothing(crs(reloaded))
                    end
                end
                # Metadata: layer attribs and per-dim attribs should match
                # after ignoring keys synthesized on write (see _LAYER_IGNORE,
                # _DIM_IGNORE). Catches real losses like dropped `bounds`
                # polygon var references.
                if exp.meta_ok
                    @test _layer_metadata_match(orig, reloaded)
                    @test _dim_metadata_match(orig, reloaded)
                else
                    @test_broken _layer_metadata_match(orig, reloaded)
                    @test_broken _dim_metadata_match(orig, reloaded)
                end
                # Every CF reference attribute (bounds, climatology,
                # coordinates, grid_mapping, geometry, cell_measures,
                # formula_terms, ancillary_variables, ...) must point at a
                # variable that actually exists in the written file. Catches
                # silent gaps where the writer preserves the attribute string
                # but never writes the referenced variable.
                dangling = _dangling_refs(outpath)
                if exp.refs_ok
                    @test isempty(dangling)
                else
                    @test_broken isempty(dangling)
                end
            end
        end
    end
end
