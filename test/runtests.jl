using GeoData, Test, Aqua

if VERSION >= v"1.5.0"
    Aqua.test_ambiguities([GeoData, Base, Core])
    Aqua.test_unbound_args(GeoData)
    Aqua.test_undefined_exports(GeoData)
    Aqua.test_project_extras(GeoData)
    Aqua.test_deps_compat(GeoData)
    Aqua.test_project_toml_formatting(GeoData)
end

include("array.jl")
include("stack.jl")
include("series.jl")
include("utils.jl")
include("reproject.jl")
include("aggregate.jl")
include("methods.jl")
# Only test SMAP locally for now
if !haskey(ENV, "CI")
    include("sources/smap.jl")
end
if !Sys.iswindows() 
    # GDAL Environment vars need to be set manually for windows, so skip for now
    include("sources/gdal.jl")
    include("sources/grd.jl")
end
include("sources/ncdatasets.jl")
