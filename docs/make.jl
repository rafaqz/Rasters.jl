using Documenter, DocumenterMarkdown
using Rasters, Logging, Statistics, Dates
using Rasters.LookupArrays, Rasters.Dimensions

# ENV["GKSwstype"] = "100"

# Plots warnings are brWarn doctests. They dont warn the second time.
# Downloads also show op in doctests. So download everything first.
function flush_info_and_warnings()
    # RasterStack(AWAP, (:tmin, :tmax); date=DateTime(2001, 1, 1))
    # RasterStack(EarthEnv{HabitatHeterogeneity}, (:evenness, :range, :contrast, :correlation))
    # RasterStack(WorldClim{BioClim}, (1, 3, 5, 7, 12))
    st = RasterStack(WorldClim{Climate}; month=1);
    #plot(st)
end
#flush_info_and_warnings()


#Logging.disable_logging(Logging.Warn)


makedocs(
    modules=[Rasters],
    clean=true,
    doctest=true,
    #format   = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    sitename="Rasters.jl",
    authors="Rafael et al.",
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

deploydocs(; repo="https://github.com/rafaqz/Rasters.jl.git", push_preview=true,
    deps=Deps.pip("mkdocs", "pygments", "python-markdown-math", "mkdocs-material",
        "pymdown-extensions", "mkdocstrings", "mknotebooks",
        "pytkdocs_tweaks", "mkdocs_include_exclude_files", "jinja2", "mkdocs-video"),
    make=() -> run(`mkdocs build`), target="site", devbranch="master")
