ENV["RASTERDATASOURCES_PATH"] = "~/Data"
using Rasters, RasterDataSources
import ArchGDAL
using CairoMakie

layers = (:evenness, :range, :contrast, :correlation)
st = RasterStack(EarthEnv{HabitatHeterogeneity}, layers)
ausbounds = X(110 .. 190), Y(-55 .. 5) # Roughly cut out Australia and New Zealand!
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
    ax = Axis(fig[1,1],
        backgroundcolor=:transparent
        )
    plt = Makie.heatmap!(ax, aus[:correlation]; colormap)
    ylims!(ax, -50, 15)
    xlims!(ax, 110, 180)
    hidedecorations!(ax)
    fig
end

mkpath(joinpath(@__DIR__, "src", "assets"))
save(joinpath(@__DIR__, "src", "assets", "logo.png"), fig; px_per_unit=2)
save(joinpath(@__DIR__, "src", "assets", "favicon.png"), fig; px_per_unit=0.25)

# do tiling for background image
function tile_background()
    grid_color = RGBAf(0.55, 0.55, 0.55, 0.16)
    fig = Figure(; figure_padding= 0, size=(600, 600),
        backgroundcolor=:transparent
        )
    ax = Axis(fig[1,1], backgroundcolor=:transparent,
        xgridcolor = grid_color, ygridcolor = grid_color,)
    hidespines!(ax, :r, :b)
    hidedecorations!(ax, grid=false)
    ax.topspinecolor = grid_color
    ax.leftspinecolor = grid_color
    ax.xticks = 1:9
    ax.yticks = 1:10
    fig
end
fig = tile_background()
save(joinpath(@__DIR__, "src", "assets", "rect_pattern.png"), fig; px_per_unit=2)
