using Rasters, RasterDataSources
import ArchGDAL
using CairoMakie

layers = (:evenness, :range, :contrast, :correlation)
st = RasterStack(EarthEnv{HabitatHeterogeneity}, layers)
ausbounds = X(110 .. 190), Y(-55 .. 15) # Roughly cut out Australia and New Zealand!
aus = st[ausbounds...] |> Rasters.trim

# colorbar attributes
colormap = :batlow
flipaxis = false
tickalign=1
width = 13
ticksize = 13
# figure
fig = with_theme(theme_dark()) do
    fig = Figure(; size=(600, 600),
        backgroundcolor=:transparent
        )
    ax = Axis(fig[1,1], #xlabel = "lon", ylabel = "lat",
        backgroundcolor=:transparent
        )
    plt = Makie.heatmap!(ax, aus[:correlation]; colormap)
    ax.xticks = 115:5:180
    ax.yticks = -50:5:5

    hidedecorations!(ax; grid=false, ticks=false)
    Label(fig[1, 1], rich("Rasters.jl"; color=:grey55, font="JuliaMono", ),
        tellwidth=false, tellheight=false,
        halign=1, valign=1, fontsize=72,)
    fig
end

mkpath(joinpath(@__DIR__, "src", "assets"))
save(joinpath(@__DIR__, "src", "assets", "logo.png"), fig; px_per_unit=2)
save(joinpath(@__DIR__, "src", "assets", "favicon.png"), fig; px_per_unit=0.25)