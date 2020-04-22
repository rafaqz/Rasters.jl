include("array.jl")
include("stack.jl")
include("series.jl")
include("aggregate.jl")
include("methods.jl")
if !Sys.iswindows() 
    # GDAL Environment vars need to be set manually for windows, so skip for now
    include("sources/grd.jl")
    include("sources/gdal.jl")
end
include("sources/ncdatasets.jl")
# Only test SMAP locally for now
# include("smap.jl")
