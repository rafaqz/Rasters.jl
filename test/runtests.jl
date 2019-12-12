include("array.jl")
include("stack.jl")
include("series.jl")
if !Sys.iswindows() 
    # GDAL Environment vars needs to be set manually for windows, so skip for now
    include("gdal.jl")
    include("grd.jl")
end
include("smap.jl")
include("ncdatasets.jl")
