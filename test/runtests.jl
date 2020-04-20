include("array.jl")
include("stack.jl")
include("series.jl")
include("aggregate.jl")
if !Sys.iswindows() 
    # GDAL Environment vars need to be set manually for windows, so skip for now
    include("gdal.jl")
    include("grd.jl")
end
include("ncdatasets.jl")
# Only test SMAP locally for now
if !(get(ENV, "CI", false))
    include("smap.jl")
end
