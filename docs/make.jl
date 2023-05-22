using Documenter, Rasters, Plots, Logging, Statistics, Dates, RasterDataSources, ArchGDAL, NCDatasets
import Makie, CairoMakie

using Rasters.LookupArrays, Rasters.Dimensions

ENV["GKSwstype"] = "100"

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
makedocs(
    modules = [Rasters],
    sitename = "Rasters.jl",
    strict = true,
    clean = false,
)

# Enable logging to console again
Logging.disable_logging(Logging.BelowMinLevel)

deploydocs(
    repo = "github.com/rafaqz/Rasters.jl.git",
)
