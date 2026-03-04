# Support for CachedCloudSource from RasterDataSources
# Enables: Raster(ERA5, :t2m) with automatic chunk caching

# Define sourcetrait for CachedCloudSource
RA.sourcetrait(::RDS.CachedCloudSource) = Zarrsource()

# Open CachedCloudSource as ZarrDataset with caching
function _open_cached_cloud(f, source::RDS.CachedCloudSource)
    store = Zarr.CachingHTTPStore(source)
    cs = Zarr.ConsolidatedStore(store, "")
    ds = ZD.ZarrDataset(cs)
    f(ds)
end

# Raster from CachedCloudSource + layer
function RA.Raster(source::RDS.CachedCloudSource, layer::Symbol;
    crs=RA.EPSG(4326),
    kw...
)
    _open_cached_cloud(source) do ds
        varname = RDS.layername(RDS.ERA5, layer)
        RA.Raster(ds[varname]; name=layer, crs, kw...)
    end
end

# Raster from ERA5 type + layer
function RA.Raster(::Type{RDS.ERA5}, layer::Symbol; kw...)
    source = RDS.getraster(RDS.ERA5)
    RA.Raster(source, layer; kw...)
end

# RasterStack from CachedCloudSource + layers
function RA.RasterStack(source::RDS.CachedCloudSource, layers::Tuple;
    crs=RA.EPSG(4326),
    kw...
)
    _open_cached_cloud(source) do ds
        rasters = map(layers) do layer
            varname = RDS.layername(RDS.ERA5, layer)
            layer => RA.Raster(ds[varname]; name=layer, crs, kw...)
        end
        RA.RasterStack(rasters...)
    end
end
function RA.RasterStack(source::RDS.CachedCloudSource, layer::Symbol; kw...)
    RA.RasterStack(source, (layer,); kw...)
end

# RasterStack from ERA5 type + layers
function RA.RasterStack(::Type{RDS.ERA5}, layers::Tuple; kw...)
    source = RDS.getraster(RDS.ERA5)
    RA.RasterStack(source, layers; kw...)
end
function RA.RasterStack(::Type{RDS.ERA5}, layer::Symbol; kw...)
    RA.RasterStack(RDS.ERA5, (layer,); kw...)
end
function RA.RasterStack(::Type{RDS.ERA5}; kw...)
    RA.RasterStack(RDS.ERA5, RDS.layers(RDS.ERA5); kw...)
end
