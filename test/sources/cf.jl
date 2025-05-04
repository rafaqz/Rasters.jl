using NearestNeighbors
# Build all ncgen files
cfdir = joinpath(dirname(dirname(Base.pathof(Rasters))), "test", "cf")
test_paths = map(filter(endswith(".ncgen"), readdir(cfdir))) do input_name
    input_path = joinpath(cfdir, input_name)
    output_path = joinpath(cfdir, "..", "data", "cf", splitext(input_name)[1] * ".nc")
    NetCDF_jll.ncgen() do exe
        run(`$exe -k nc4 -b -o $output_path $input_path`)
    end
    output_path
end
rasters = map(test_paths) do test_path
    name = splitext(basename(test_path))[1]
    try
        name => Raster(test_path; lazy=true)
    catch
        name => nothing
    end
end 
rasters[1][2]
rasterstacks = map(test_paths) do test_path
    name = splitext(basename(test_path))[1]
    name == "5_1_independent_coordinate_variables" && return name => nothing
    name == "5_10_british_national_grid" && return name => nothing
    println(name)
    name => RasterStack(test_path; lazy=true)
end ;
read(last.(rasterstacks)[1].str_variable)