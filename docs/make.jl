using Documenter, Rasters, Plots, Logging, Statistics, Dates, 
    RasterDataSources, ArchGDAL, NCDatasets, CoordinateTransformations
import Makie, CairoMakie
using DocumenterVitepress
using Rasters.LookupArrays, Rasters.Dimensions
import Shapefile, DataFrames, NaturalEarth # to avoid precompilation in doctests

# Don't output huge svgs for Makie plots
CairoMakie.activate!(type = "png")

# Plots warnings are brWarn doctests. They don't warn the second time.
# Downloads also show op in doctests. So download everything first.
function flush_info_and_warnings()
    RasterStack(AWAP, (:tmin, :tmax); date=DateTime(2001, 1, 1))
    RasterStack(EarthEnv{HabitatHeterogeneity}, (:evenness, :range, :contrast, :correlation))
    RasterStack(WorldClim{BioClim}, (1, 3, 5, 7, 12))
    st = RasterStack(WorldClim{Climate}; month=1);
    plot(st)
end
flush_info_and_warnings()


Logging.disable_logging(Logging.Warn)

# Make the docs, without running the tests again
# We need to explicitly add all the extensions here

makedocs(
    modules = [
        Rasters,
        Base.get_extension(Rasters, :RastersArchGDALExt),
        Base.get_extension(Rasters, :RastersCoordinateTransformationsExt),
        Base.get_extension(Rasters, :RastersMakieExt),
        Base.get_extension(Rasters, :RastersNCDatasetsExt),
        Base.get_extension(Rasters, :RastersRasterDataSourcesExt),
    ],
    sitename = "Rasters.jl",
    authors="Rafael Schouten et al.",
    clean=true,
    doctest=true,
    checkdocs=:all,
    format=DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/rafaqz/Rasters.jl", # this must be the full URL!
        devbranch = "main",
        devurl = "dev";
    ),
    source = "src",
    build = "build",
    warnonly=false,
)

# Enable logging to console again
Logging.disable_logging(Logging.BelowMinLevel)

DocumenterVitepress.deploydocs(; repo="github.com/rafaqz/Rasters.jl",
    branch = "gh-pages",
    devbranch = "main",
    push_preview = true
    )
