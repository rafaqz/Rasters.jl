using Rasters, NCDatasets 
using NCDatasets.NetCDF_jll
using NearestNeighbors
using OrderedCollections
using Rasters.Lookups
using Test

# Build all ncgen files
cfdir = joinpath(dirname(dirname(Base.pathof(Rasters))), "test", "cf")
example_paths = map(filter(endswith(".ncgen"), readdir(cfdir))) do input_name
    input_path = joinpath(cfdir, input_name)
    output_path = joinpath(cfdir, "..", "data", "cf", splitext(input_name)[1] * ".nc")
    NetCDF_jll.ncgen() do exe
        run(`$exe -k nc4 -b -o $output_path $input_path`)
    end
    output_path
end
# rasters = map(example_paths) do example_path
#     name = splitext(basename(example_path))[1]
#     name == "5_1_independent_coordinate_variables" && return name => nothing
#     name => Raster(example_path; lazy=true)
# end 
# rasterstacks = map(enumerate(example_paths)) do (i, example_path)
#     name = splitext(basename(example_path))[1]
#     println(i, " => ", name)
#     name == "5_1_independent_coordinate_variables" && return name => nothing
#     name => RasterStack(example_path; lazy=true)
# end ;
# st = last.(rasterstacks)[12]
# st = last.(rasterstacks)[14]
# refdims(st)
# inds = findall(.!(isempty.(refdims.(last.(rasterstacks)))))
# refdims(last.(rasterstacks[inds])[1])

# for (k, r) in rasters
#     println(k)
#     display(r)
#     readline()
# end

example_ids = replace.(strip.(first.(split.(basename.(example_paths), r"[A-Za-z]")), '_'), ("_" => ".",))
examples = OrderedDict(example_ids .=> example_paths)

@testset "2.1" begin
    rast = RasterStack(examples["2.1"])
    @test rast.char_variable == rast.str_variable
    # Also works lazily
    @test read(RasterStack(examples["2.1"]; lazy=true).char_variable) == 
          read(RasterStack(examples["2.1"]; lazy=true)).char_variable == 
          rast.str_variable
    # Fails due to DiskArrays.jl size checks during the comparison iteration
    @test_broken RasterStack(examples["2.1"]; lazy=true).char_variable == rast.str_variable
end

@testset "3.1" begin
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

@testset "5.1" begin
    rast = RasterStack(examples["5.1"])
    @test dims(rast) == (X(NoLookup(1:36)), Y(NoLookup(1:18)), Dim{:pres}(NoLookup(1:15)), Ti(NoLookup(1:4)))
    @test metadata(rast.xwind)["long_name"] == "zonal wind"
    @test metadata(rast.xwind)["units"] == "m/s"
    # TODO add real vars to the file so there is also dimension metadata
end

rast = RasterStack(examples["5.2"])
rast = RasterStack(examples["5.3"])
rast = RasterStack(examples["5.6"])
rast = RasterStack(examples["5.8"])

