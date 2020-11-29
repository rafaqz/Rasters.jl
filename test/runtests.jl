using GeoData, Test, Aqua, SafeTestsets

if VERSION >= v"1.5.0"
    #Aqua.test_ambiguities([GeoData, Base, Core])
    Aqua.test_unbound_args(DimensionalData)
    Aqua.test_undefined_exports(DimensionalData)
    Aqua.test_project_extras(DimensionalData)
    Aqua.test_stale_deps(DimensionalData)
    Aqua.test_deps_compat(DimensionalData)
    Aqua.test_project_toml_formatting(DimensionalData)
    Aqua.test_project_extras(DimensionalData)
    Aqua.test_stale_deps(DimensionalData)
end

@time @safetestset "array" begin include("array.jl") end
@time @safetestset "stack" begin include("stack.jl") end
@time @safetestset "series" begin include("series.jl") end
@time @safetestset "utils" begin include("utils.jl") end
@time @safetestset "set" begin include("set.jl") end
@time @safetestset "reproject" begin include("reproject.jl") end
@time @safetestset "aggregat" begin include("aggregate.jl") end
@time @safetestset "methods" begin include("methods.jl") end
@time @safetestset "resample" begin include("resample.jl") end
# Only test SMAP locally for now
if !haskey(ENV, "CI")
    @time @safetestset "smap" begin include("sources/smap.jl") end
end
if !Sys.iswindows()
    # GDAL Environment vars need to be set manually for windows, so skip for now
    @time @safetestset "gdal" begin include("sources/gdal.jl") end
    @time @safetestset "grd" begin include("sources/grd.jl") end
end
@time @safetestset "ncdatasets" begin include("sources/ncdatasets.jl") end
