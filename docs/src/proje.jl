ENV["RASTERDATASOURCES_PATH"] = "/Users/lalonso/Data"
using Rasters, RasterDataSources, ArchGDAL
using DimensionalData
using DimensionalData.Lookups
using NaNStatistics
using GLMakie

# resample section
ras = Raster(WorldClim{BioClim}, 5)
ras_m = replace_missing(ras, missingval=NaN) # ArchGDAL warp handles NaNs via the metadata in the Raster

ras_sample = resample(ras_m; size=(2160, 1080))

methods = ["average", "mode", "max", "sum"]
sizes = [(2160, 1080), (1440, 720), (720, 360), (360, 180)]
resolutions = [0.16666666666666666, 0.25, 0.5, 1.0]

# other available methods to try:  "bilinear", "cubic", "cubicspline", "lanczos", "min", "med", "q1", "q3", "near".

method_sizes = [resample(ras_m; size=size, method=method) for method in methods for size in sizes]
method_res = [resample(ras_m; res=res, method=method) for method in methods for res in resolutions]

with_theme(Rasters.theme_rasters()) do
    colorrange = (nanminimum(ras_m), nanmaximum(ras_m))
    hm=nothing
    fig = Figure(; size = (1000, 600))
    axs = [Axis(fig[i,j], title="size=$(size), method=:$(method)", titlefont=:regular)
        for (i, method) in enumerate(methods) for (j, size) in enumerate(sizes)]
    for (i, ax) in enumerate(axs)
        hm = heatmap!(ax, method_sizes[i]; colorrange)
    end
    Colorbar(fig[:,end+1], hm)
    hidedecorations!.(axs; grid=false)
    rowgap!(fig.layout, 5)
    colgap!(fig.layout, 10)
    fig
end

with_theme(Rasters.theme_rasters()) do
    colorrange = (nanminimum(ras_m), nanmaximum(ras_m))
    hm=nothing
    fig = Figure(; size = (1000, 600))
    axs = [Axis(fig[i,j], title="res=$(round(res, digits=4)), method=:$(method)", titlefont=:regular)
        for (i, method) in enumerate(methods) for (j, res) in enumerate(resolutions)]
    for (i, ax) in enumerate(axs)
        hm = heatmap!(ax, method_res[i]; colorrange)
    end
    Colorbar(fig[:,end+1], hm)
    hidedecorations!.(axs; grid=false)
    rowgap!(fig.layout, 5)
    colgap!(fig.layout, 10)
    fig
end

# SINUSOIDAL_CRS = ProjString("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

SINUSOIDAL_CRS = ProjString("+proj=sinu +lon_0=0 +type=crs")

fig, ax, plt = heatmap(ras_m)
Colorbar(fig[1,2], plt)
fig

ras_sin = resample(ras_m; size=(2160, 1080), crs=SINUSOIDAL_CRS, method="average")

fig, ax, plt = heatmap(ras_sin)
Colorbar(fig[1,2], plt)
display(GLMakie.Screen(), fig)


ras_epsg = resample(ras_sin; size=(1440,720), crs=EPSG(4326), method="average")

fig = heatmap(ras_epsg)
display(GLMakie.Screen(), fig)

locus_resampled = DimensionalData.shiftlocus(Center(), ras_epsg)
heatmap(locus_resampled)

x_range = LinRange(-180, 179.75, 1440)
y_range = LinRange(89.75, -90, 720)
ras_data = ras_epsg.data

ras_scratch = Raster(ras_data, (X(x_range; sampling=Intervals(Start())),
    Y(y_range; sampling=Intervals(Start()))), crs=EPSG(4326), missingval=NaN)
heatmap(ras_scratch)
# locus_resampled = DimensionalData.shiftlocus(Center(), ras_scratch)

# SINUSOIDAL_CRS = ProjString("+proj=sinu +R=6371007.181 +nadgrids=@null +wktext") # here, both output same bounds!
ras_sin_s = resample(ras_scratch; size=(1440,720), crs=SINUSOIDAL_CRS, method="average")
heatmap(ras_sin_s)

ras_sin_2 = resample(ras_epsg; size=(1440,720), crs=SINUSOIDAL_CRS, method="average")
heatmap(ras_sin_2)

ras_epsg = resample(ras_sin_s; size=(1440,720), crs=EPSG(4326), method="average")
locus_resampled = DimensionalData.shiftlocus(Center(), ras_epsg)

heatmap(locus_resampled)