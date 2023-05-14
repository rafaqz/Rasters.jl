using Rasters, Test, Aqua, SafeTestsets

if VERSION >= v"1.5.0"
    # Aqua.test_ambiguities([Rasters, Base, Core])
    Aqua.test_unbound_args(Rasters)
    Aqua.test_stale_deps(Rasters)
    Aqua.test_undefined_exports(Rasters)
    Aqua.test_project_extras(Rasters)
    Aqua.test_deps_compat(Rasters)
    Aqua.test_project_toml_formatting(Rasters)
end

@time @safetestset "array" begin include("array.jl") end
@time @safetestset "stack" begin include("stack.jl") end
@time @safetestset "series" begin include("series.jl") end
@time @safetestset "utils" begin include("utils.jl") end
@time @safetestset "set" begin include("set.jl") end
@time @safetestset "reproject" begin include("reproject.jl") end
@time @safetestset "aggregate" begin include("aggregate.jl") end
@time @safetestset "rasterize" begin include("rasterize.jl") end
@time @safetestset "methods" begin include("methods.jl") end
@time @safetestset "adapt" begin include("adapt.jl") end
# Only test SMAP locally for now, also RasterDataSources because CI dowloads keep breaking
if !haskey(ENV, "CI")
    @time @safetestset "smap" begin include("sources/smap.jl") end
    @time @safetestset "rasterdatasources" begin include("sources/rasterdatasources.jl") end
end
if !Sys.iswindows()
    # GDAL Environment vars need to be set manually for windows, so skip for now
    @time @safetestset "gdal" begin include("sources/gdal.jl") end
    @time @safetestset "grd" begin include("sources/grd.jl") end
end
@time @safetestset "ncdatasets" begin include("sources/ncdatasets.jl") end
@time @safetestset "plot recipes" begin include("plotrecipes.jl") end
