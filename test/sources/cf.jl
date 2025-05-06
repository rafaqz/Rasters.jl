using Rasters, NCDatasets 
using NCDatasets.NetCDF_jll
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
    name == "5_1_independent_coordinate_variables" && return name => nothing
    name => Raster(test_path; lazy=true)
end 
rasterstacks = map(enumerate(test_paths)) do (i, test_path)
    name = splitext(basename(test_path))[1]
    println(i, " => ", name)
    name == "5_1_independent_coordinate_variables" && return name => nothing
    name => RasterStack(test_path; lazy=true)
end ;
st = last.(rasterstacks)[12]
st = last.(rasterstacks)[14]
refdims(st)
inds = findall(.!(isempty.(refdims.(last.(rasterstacks)))))
refdims(last.(rasterstacks[inds])[1])

for (k, r) in rasters
    println(k)
    display(r)
    readline()
end