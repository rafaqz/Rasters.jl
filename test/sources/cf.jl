using Rasters, NCDatasets 
using NCDatasets.NetCDF_jll
using NearestNeighbors
using OrderedCollections
using Rasters.Lookups
using Test
using Rasters: name
using Dates
using GeoInterface

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

example_ids = replace.(strip.(first.(split.(basename.(example_paths), r"[A-Za-z]")), '_'), ("_" => ".",))
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
    # Also works lazily
    @test read(RasterStack(examples["2.1"]; lazy=true).char_variable) == 
          read(RasterStack(examples["2.1"]; lazy=true)).char_variable == 
          rast.str_variable
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
    rast = RasterStack(examples["5.1"])
    @test dims(rast) == (X(NoLookup(1:36)), Y(NoLookup(1:18)), Dim{:pres}(NoLookup(1:15)), Ti(NoLookup(1:4)))
    @test metadata(rast.xwind)["long_name"] == "zonal wind"
    @test metadata(rast.xwind)["units"] == "m/s"
    # TODO add real vars to the file so there is also dimension metadata
end

@testset "5.2 two-dimensional coordinate variables" begin
    rast = RasterStack(examples["5.2"]; lazy=true)
    # TODO multi-dimensional mapped
    @test_broken rast[Lat=1.0, Lon=2.0] isa Float32
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
    @test_broken rast[Lat=1.0, Lon=2.0] isa Float32
end

@testset "5.8-9 WGS 84 EPSG" begin
    rast = RasterStack(examples["5.8"]; lazy=true)
    @test crs(rast) isa EPSG
    @test crs(rast) == EPSG(4326)
    rast = RasterStack(examples["5.9"]; lazy=true)
    @test crs(rast) isa EPSG
    @test crs(rast) == EPSG(4326)
end

@testset "5.10 British national grid" begin
    rast = RasterStack(examples["5.10"]; lazy=true)
    @test_broken !isnothing(crs(rast))
    @test keys(rast) == (:temp, :pres)
end

@testset "5.11 WGS 84 WellKnownText2" begin
    rast = RasterStack(examples["5.11"]; lazy=true)
    @test crs(rast) isa WellKnownText2
end

@testset "5.14 refdims from scalar coordinates" begin
    rast = RasterStack(examples["5.14"]; lazy=true)
    @test name(refdims(rast)) == (:atime, :p500)
end

@testset "domains" begin
    RasterStack(examples["5.15"]; lazy=true)
    RasterStack(examples["5.16"]; lazy=true)
    @testset "5.17 domain defines dimension and coords" begin
        RasterStack(examples["5.17"]; lazy=true)
    end
    RasterStack(examples["5.18"]; lazy=true)
    RasterStack(examples["5.19"]; lazy=true)
    RasterStack(examples["5.20"]; lazy=true)
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
    RasterStack(examples["7.2"]; lazy=true)
end

@testset "7.3-4 formula terms" begin
    # These are not implemented
    RasterStack(examples["7.3"]; lazy=true)
    RasterStack(examples["7.4"]; lazy=true)
end

@testset "7.5 methods applied to a timeseries" begin
    rast = RasterStack(examples["7.5"]; lazy=true)
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
    # Not implemented
    RasterStack(examples["7.9"]; lazy=true)
    RasterStack(examples["7.10"]; lazy=true)
    RasterStack(examples["7.11"]; lazy=true)
    RasterStack(examples["7.12"]; lazy=true)
    RasterStack(examples["7.13"]; lazy=true)
    RasterStack(examples["7.14"]; lazy=true)
end

@testset "Geometry lookups" begin
    @testset "7.15 linestring" begin
        rast = RasterStack(examples["7.15"])
        @test size(rast) == (4, 2)
        @test dims(rast, :Geometry)[1] isa GeoInterface.Wrappers.LineString
        @test GeoInterface.getpoint(dims(rast, :Geometry)[1]) == [(30.0, 10.0), (10.0, 30.0), (40.0, 40.0)]
    end
    @testset "7.16 polygon" begin
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
