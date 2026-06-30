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
    rast = RasterStack(examples["2.1"])
    expected = Union{Missing,String}["alpha", "beta", "gamma", "delta", "epsilon", "zeta"]
    @test parent(rast.char_variable) == expected
    @test parent(rast.str_variable) == expected
    @test rast.char_variable == rast.str_variable
    @test metadata(rast.char_variable)["long_name"] == "strings of type char"
    @test metadata(rast.str_variable)["long_name"] == "strings of type string"
    # Round-trip through write/read preserves equality of both representations
    mktempdir() do dir
        path = joinpath(dir, "string_roundtrip.nc")
        write(path, rast; force=true)
        rast2 = RasterStack(path)
        @test parent(rast2.char_variable) == expected
        @test parent(rast2.str_variable) == expected
        @test rast2.char_variable == rast2.str_variable
    end
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
    @test keys(rast) == (:T,)
    @test size(rast) == (4, 3, 2)
    @test name.(dims(rast)) == (:X, :Y, :Z)
    # 2D auxiliary lon/lat coords produce ProjectedArrayLookup on X and Y
    @test lookup(rast, X) isa Rasters.ProjectedArrayLookup
    @test lookup(rast, Y) isa Rasters.ProjectedArrayLookup
    @test size(lookup(rast, X).matrix) == (4, 3)
    @test size(lookup(rast, Y).matrix) == (4, 3)
    # Z is a regular 1D pressure axis
    @test lookup(rast, Z) isa Sampled
    @test collect(lookup(rast, Z)) == Float32[1000.0, 500.0]
    @test metadata(dims(rast, Z))["units"] == "hPa"
    @test metadata(dims(rast, Z))["long_name"] == "pressure level"
    @test metadata(rast.T)["units"] == "K"
    @test metadata(rast.T)["long_name"] == "temperature"
    # Integer indexing through 3D data: T values are 1..24 in order
    @test rast.T[X=1, Y=1, Z=1] === 1.0f0
    @test rast.T[X=4, Y=3, Z=2] === 24.0f0
    # Near() via the 2D lon/lat coordinate matrix
    # lon[X=1,Y=1]=10, lat[X=1,Y=1]=50 -> T=1
    @test rast.T[X=Near(10.0), Y=Near(50.0), Z=1] === 1.0f0
    # lon[X=4,Y=1]=40, lat[X=4,Y=1]=53 -> T=4
    @test rast.T[X=Near(40.0), Y=Near(53.0), Z=1] === 4.0f0
end

@testset "5.3 reduced horizontal grid" begin
    rast = RasterStack(examples["5.3"]; lazy=true)
    @test rast[rgrid=2].PS === 0.2f0
    @test rast[rgrid=At(3.0, 30.0)].PS === 0.3f0
    @test rast[X=At(4.0)].PS == Union{Missing,Float32}[0.4f0]
    @test rast[Y=At(30.0)].PS == Union{Missing,Float32}[0.3f0]
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
    # CRS is now parsed via CFCoordinateReferenceSystems
    # This file has multiple CRS mappings (crsOSGB for x/y, crsWGS84 for lat/lon)
    # so it returns a Vector of Pairs
    c = crs(rast)
    @test !isnothing(c)
    @test c isa Vector{<:Pair}
    @test first(c).first isa WellKnownText2
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
        # TODO : cell_measures (cell_area attached to a domain) is not yet a first-class
        # concept in Rasters. The domain variable is dropped, cell_area is loaded as
        # a regular layer, and the cell dim merges lon/lat into a MergedLookup.
        rast = RasterStack(examples["5.17"]; lazy=true)
        @test keys(rast) == (:cell_area,)
        @test name.(dims(rast)) == (:cell, :Ti)
        @test lookup(rast, :cell) isa MergedLookup
        @test name.(lookup(rast, :cell).dims) == (:X, :Y)
        @test length(lookup(rast, :cell)) == 2562
        @test length(lookup(rast, Ti)) == 12
        @test metadata(rast.cell_area)["units"] == "m2"
        @test metadata(rast.cell_area)["standard_name"] == "cell_area"
        @test metadata(rast.cell_area)["long_name"] == "area of grid cell"
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

@testset "5.20 indexed ragged array" begin
    # CF indexed ragged array (Chapter 9): stationIndex(obs) maps each obs to a
    # station, so the (obs,) flat axis really encodes (station, time) pairs.
    # Current Rasters: stack is loaded flat - humidity(obs), station_info(station)
    # are separate layers, no time lookup is built from stationIndex+time.
    # Ideal: humidity should be unpacked to a (station, time) Raster, with
    # station coords (lon, lat, alt, station_name) and a Ti lookup.
    rast = RasterStack(examples["5.20"]; lazy=true)
    @test :humidity in keys(rast)
    @test :station_info in keys(rast)
    @test parent(rast.station_info) == [1, 2, 3]
    @test metadata(rast.humidity)["standard_name"] == "relative_humidity"
    @test metadata(rast.humidity)["units"] == "%"
    @test name.(dims(rast)) == (:station, :obs)
    # --- broken: ragged-array unpacking is not implemented ---
    # humidity should be unpacked into a 2D (station, time) raster
    @test_broken size(rast.humidity) == (3, 4)
    @test_broken hasdim(rast.humidity, Ti)
    @test_broken hasdim(rast.humidity, :station)
    # Time values should be parsed into a Ti lookup
    @test_broken lookup(rast, Ti) == DateTime(1970, 1, 1):Day(1):DateTime(1970, 1, 4)
    # Per-station coords should attach to the station dim
    @test_broken lookup(rast, :station) isa MergedLookup
    # The stationIndex helper variable should not be exposed as a layer once unpacked
    @test_broken !(:stationIndex in keys(rast))
    # humidity values for station 1 (Station_A) should be [50, 51, 52, 53]
    @test_broken parent(rast.humidity[station=1]) == Float32[50, 51, 52, 53]
end

@testset "5.21 mesh" begin
    # UGRID mesh topology is not yet interpreted by Rasters - the stack loads as
    # a flat collection of layers/dims. Skipping deeper assertions until a
    # MeshLookup (or similar) representation is designed.
    rast = RasterStack(examples["5.21"]; lazy=true)
    @test :mesh in keys(rast)
    @test metadata(rast.mesh)["cf_role"] == "mesh_topology"
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
    # Selection by direct inner dimension name
    @test size(rast[Dim{:taxon_lsid}(At("urn:lsid:marinespecies.org:taxname:104464"))]) == (1, 100)
    @test size(rast[Dim{:taxon_name}(At("Calanus finmarchicus"))]) == (1, 100)
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
    @test size(rast[Dim{:model_level}(At(Int32(10)))]) == (30, 1)
    @test size(rast[Dim{:model_level}(At(Int32(30)))]) == (30, 1)
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

@testset "7.4 cell areas" begin
    # PS has cell_measures = "area: cell_area". Currently cell_area is loaded as
    # a sibling layer. Long-term it should be attached to PS (cell_measures support).
    rast = RasterStack(examples["7.4"]; lazy=true)
    @test keys(rast) == (:PS, :cell_area)
    @test size(rast.PS) == (2562, 12)
    @test size(rast.cell_area) == (2562,)
    @test name.(dims(rast)) == (:cell, :Ti)
    @test lookup(rast, :cell) isa MergedLookup
    @test length(lookup(rast, :cell)) == 2562
    @test length(lookup(rast, Ti)) == 12
    @test metadata(rast.PS)["units"] == "Pa"
    @test metadata(rast.PS)["cell_measures"] == "area: cell_area"
    @test metadata(rast.cell_area)["units"] == "m2"
    @test metadata(rast.cell_area)["standard_name"] == "cell_area"
    # cell_measures should ultimately attach cell_area to PS rather than expose it as a sibling layer
    @test_broken !(:cell_area in keys(rast))
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
