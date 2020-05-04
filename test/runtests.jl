include("array.jl")
include("stack.jl")
include("series.jl")
include("aggregate.jl")
include("methods.jl")
# Only test SMAP locally for now
if !haskey(ENV, "CI") || !ENV["CI"] 
    include("sources/smap.jl")
end
if !Sys.iswindows() 
    # GDAL Environment vars need to be set manually for windows, so skip for now
    include("sources/gdal.jl")
    include("sources/grd.jl")
end
include("sources/ncdatasets.jl")
