using Rasters, RasterDataSources, ArchGDAL
using DimensionalData
using DimensionalData.Lookups
using NaNStatistics
using GLMakie

ras = Raster(WorldClim{BioClim}, 5)
ras_m = replace_missing(ras, missingval=NaN)
SINUSOIDAL_CRS = ProjString("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs")

ras_sin = resample(ras_m; size=(2160, 1080), crs=SINUSOIDAL_CRS, method="average")
heatmap(ras_sin)

ras_epsg = resample(ras_sin; size=(1440,720), crs=EPSG(4326), method="average")
locus_resampled = DimensionalData.shiftlocus(Center(), ras_epsg)


x_range = LinRange(-180, 179.75, 1440)
y_range = LinRange(89.75, -90, 720)
ras_data = ras_025.data

ras_scratch = Raster(ras_data, (X(x_range; sampling=Intervals(Start())),
    Y(y_range; sampling=Intervals(Start()))), crs=EPSG(4326))
heatmap(ras_scratch)
# locus_resampled = DimensionalData.shiftlocus(Center(), ras_scratch)

# SINUSOIDAL_CRS = ProjString("+proj=sinu +R=6371007.181 +nadgrids=@null +wktext") # here, both output same bounds!
ras_sin_s = resample(ras_scratch; size=(1440,720), crs=SINUSOIDAL_CRS, method="average")
heatmap(ras_sin_s)

ras_sin_2 = resample(ras_025; size=(1440,720), crs=SINUSOIDAL_CRS, method="average")
heatmap(ras_sin_2)

ras_epsg = resample(ras_sin_s; size=(1440,720), crs=EPSG(4326), method="average")
locus_resampled = DimensionalData.shiftlocus(Center(), ras_epsg)

heatmap(locus_resampled)