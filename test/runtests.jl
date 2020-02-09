if !Sys.iswindows() 
    # GDAL Environment vars need to be set manually for windows, so skip for now
    include("grd.jl")
    include("gdal.jl")
end
include("array.jl")
include("stack.jl")
include("series.jl")
include("smap.jl")
include("ncdatasets.jl")
