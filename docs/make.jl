using Documenter, Rasters, Plots, Logging, Statistics, Dates, 
    RasterDataSources, ArchGDAL, NCDatasets, HDF5, CoordinateTransformations
import Makie, CairoMakie
using DocumenterMarkdown
using Rasters.LookupArrays, Rasters.Dimensions

ENV["GKSwstype"] = "100"
ENV["RASTERDATASOURCES_PATH"] = "/Users/lalonso/Data/"

# Plots warnings are brWarn doctests. They dont warn the second time.
# Downloads also show op in doctests. So download everything first.
function flush_info_and_warnings()
    # RasterStack(AWAP, (:tmin, :tmax); date=DateTime(2001, 1, 1))
    # RasterStack(EarthEnv{HabitatHeterogeneity}, (:evenness, :range, :contrast, :correlation))
    # RasterStack(WorldClim{BioClim}, (1, 3, 5, 7, 12))
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
        Base.get_extension(Rasters, :RastersHDF5Ext),
        Base.get_extension(Rasters, :RastersMakieExt),
        Base.get_extension(Rasters, :RastersNCDatasetsExt),
        Base.get_extension(Rasters, :RastersRasterDataSourcesExt),
    ],
    sitename = "Rasters.jl",
    authors="Rafael Schouten et al.",
    clean=true,
    doctest=true,
    strict=[
        :doctest,
        :linkcheck,
        :parse_error,
        :example_block,
        # Other available options are
        # :autodocs_block, :cross_references, :docs_block, :eval_block, :example_block,
        # :footnote, :meta_block, :missing_docs, :setup_block
    ], checkdocs=:all, format=Markdown(), draft=false,
    build=joinpath(@__DIR__, "docs")
)

# Enable logging to console again
Logging.disable_logging(Logging.BelowMinLevel)

deploydocs(; repo="github.com/rafaqz/Rasters.jl.git", push_preview=true,
    deps=Deps.pip("mkdocs", "pygments", "python-markdown-math", "mkdocs-material",
        "pymdown-extensions", "mkdocstrings", "mknotebooks",
        "pytkdocs_tweaks", "mkdocs_include_exclude_files", "jinja2", "mkdocs-video"),
    make=() -> run(`mkdocs build`),
    versions = ["stable" => "v^", "v#.#.#", "dev" => "dev"],
    target="site")
