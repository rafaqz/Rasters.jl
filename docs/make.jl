using Documenter, GeoData, Plots, Logging

ENV["GKSwstype"] = "100"

# Disable logging as @info messes with doctests
Logging.disable_logging(Logging.Warn)

# Make the docs, without running the tests again
makedocs(
    modules = [GeoData],
    sitename = "GeoData.jl",
    strict = false,
)

# Enable logging to console again
Logging.disable_logging(Logging.BelowMinLevel)

deploydocs(
    repo = "github.com/rafaqz/GeoData.jl.git",
)
