include("array.jl")
include("stack.jl")
include("series.jl")
# GDAL Environment vars needs to be set manually for windows, so skip for now
Sys.iswindows() || include("gdal.jl")
include("ncdatasets.jl")
include("smap.jl")
include("grd.jl")
