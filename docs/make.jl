using Documenter, GeoData, Plots, Logging, Statistics, Dates

using GeoData.LookupArrays, GeoData.Dimensions

ENV["GKSwstype"] = "100"

# Plots warnings are brWarn doctests. They dont warn the second time.
# Downloads also show op in doctests. So download everything first.
function flush_info_and_warnings()
    # GeoStack(AWAP, (:tmin, :tmax); date=DateTime(2001, 1, 1))
    # GeoStack(EarthEnv{HabitatHeterogeneity}, (:evenness, :range, :contrast, :correlation))
    # GeoStack(WorldClim{BioClim}, (1, 3, 5, 7, 12))
    st = GeoStack(WorldClim{Climate}; month=1);
    plot(st)
end
flush_info_and_warnings()


Logging.disable_logging(Logging.Warn)

# Make the docs, without running the tests again
makedocs(
    modules = [GeoData],
    sitename = "GeoData.jl",
    strict = true,
    clean = false,
)

# Enable logging to console again
Logging.disable_logging(Logging.BelowMinLevel)

deploydocs(
    repo = "github.com/rafaqz/GeoData.jl.git",
)
