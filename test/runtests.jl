include("array.jl")
include("stack.jl")
include("series.jl")
include("ncdatasets.jl")
include("smap.jl")
if !Sys.iswindows() 
    # GDAL Environment vars needs to be set manually for windows, so skip for now
    include("gdal.jl")
    include("grd.jl")
end
