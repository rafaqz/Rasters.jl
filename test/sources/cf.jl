using Rasters
import NCDatasets
import NCDatasets.NetCDF_jll
using NearestNeighbors
using OrderedCollections
using Rasters.Lookups
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
    @test rast.char_variable == rast.str_variable
    # Fails due to DiskArrays.jl size checks during the comparison iteration
    @test_broken RasterStack(examples["2.1"]; lazy=true).char_variable == rast.str_variable
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
    @test rast[X=Near(1.0), Y=Near(2.0), Z=1].T isa Float32
end

@testset "5.3 reduced horizontal grid" begin
    rast = RasterStack(examples["5.3"]; lazy=true)
    @test rast[rgrid=2].PS === 0.2f0
    @test rast[rgrid=At(3.0, 30.0)].PS === 0.3f0
    @test_broken rast[X=At(4.0)].PS == [0.4]
    @test_broken rast[Y=At(30.0)].PS == [0.3]
    @test_broken rast[X=At(1.0), Y=At(10.0)] == 0.1
end

@testset "5.6 rotated pole grid" begin
    rast = RasterStack(examples["5.6"]; lazy=true)
    @test lookup(rast, :rlat) isa ArrayLookup
    @test lookup(rast, :rlon) isa ArrayLookup
    @test lookup(rast, Z) == DimensionalData.format(Sampled([100.0f0, 200.0f0]; span=Regular(100.0f0)))
    @test_broken rast.temp[X=Near(-0.9), Y=Near(3.0), Z=1] === 3.0f0
    @test_broken rast.temp[X=At(-0.44758424f0), Y=At(2.1773858f0), Z=1] === 2.0f0
end

@testset "5.8-9 WGS 84 EPSG" begin
    # TODO this is wrong
    rast = RasterStack(examples["5.8"]; lazy=true)
    @test crs(rast) isa EPSG
    @test crs(rast) == EPSG(4326)
    rast = RasterStack(examples["5.9"]; lazy=true)
    @test crs(rast) isa EPSG
    @test crs(rast) == EPSG(4326)
end

@testset "5.10 British national grid" begin
    rast = RasterStack(examples["5.10"]; lazy=true);
    @test keys(rast) == (:temp, :pres)
    # Need CFCRS.jl for this
    @test_broken !isnothing(crs(rast))
end

@testset "5.11 WGS 84 WellKnownText2" begin
    rast = RasterStack(examples["5.11"]; lazy=true)
    # TODO extent is broken
    # @test extent(rast)
    @test crs(rast) isa WellKnownText2
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
        @test dims(rast) isa Tuple{<:Z,<:Dim{:rlat},<:Dim{:rlon}}
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

@testset "taxon name/id" begin
    rast = Raster(examples["6.1.2"]; lazy=true)
    @test name(rast) == :abundance
    @test lookup(rast, :taxon_lsid) isa Categorical{String,Vector{String},Unordered}
    @test lookup(rast, :taxon_lsid)[1] == "urn:lsid:marinespecies.org:taxname:104464"
    @test metadata(lookup(rast, :taxon_lsid))["standard_name"] == "biological_taxon_lsid"
    # Second lookup not yet implemented
    @test_broken lookup(rast, :taxon_name)[1] == "Calanus finmarchicus"
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

@testset "6.2 seondary lookup coordinates" begin
    rast = Raster(examples["6.2"]; lazy=true)
    @test name(rast) == :xwind
    @test lookup(rast, :sigma) == -1:-1:-5
    # Second lookup not yet implemented
    @test_broken lookup(rast, :model_level) == 10:10:50
end

@testset "7.2 non-aligned horizontal grid" begin
    # TODO: load bounds as squarish polygons
    rast = RasterStack(examples["7.2"]; lazy=true)
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
    @test rast[Ti=Near(DateTime(1960, 2))] == rast[Ti=1]
    @test rast[Ti=At(DateTime(1960, 4, 16))] ==
          rast[Ti=At(DateTime(1970, 4, 16))] ==
          rast[Ti=At(DateTime(1990, 4, 16))] == rast[Ti=1]
    @test_throws Lookups.SelectorError rast[Ti=At(DateTime(2000, 4, 16))]
    @test_throws Lookups.SelectorError rast[Ti=At(DateTime(1950, 4, 16))]
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
