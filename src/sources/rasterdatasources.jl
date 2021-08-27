using .RasterDataSources

using .RasterDataSources: RasterDataSource

const RDS = RasterDataSources

const LayerItr = Union{AbstractArray,Tuple}

"""
    GeoArray(T::Type{<:RasterDataSource}, [layer]; kw...) => AbstractArray

Load a `RasterDataSource` as an `AbstractGeoArray`. `T`, `args` are
are passed to `getraster`, while `kw` args are for both `getraster` and
`AbstractGeoArray`.

# Keywords

- `month`: For `Climate` datasets
- `date`: For `Weather` datasets
- `res`: For datasets with multiple resolutions

Normal `GeoArray` keywords are passed to the constructor.
"""
GeoArray(T::Type{<:RasterDataSource}; kw...) = GeoArray(T, first(RDS.layers(T)); kw...) 
function GeoArray(T::Type{<:RasterDataSource}, layer; kw...)
    rds_kw, gd_kw = _filterkw(kw)
    filename = getraster(T, layer; rds_kw...)
    GeoArray(filename; name=_layerkey(T, layer), _sourcekw(T)..., gd_kw...)
end

"""
    GeoStack(T::Type{<:RasterDataSource}, [layers::Union{Symbol,AbstractArray,Tuple}]; kw...) => AbstractGeoStack

Load a `RasterDataSource` as an `AbstractGeoStack`. `T`, `args` are
are passed to `getraster`, while `kw` args are for both `getraster` and
`AbstractGeoStack`.

# Keywords

- `month`: For `Climate` datasets
- `date`: For `Weather` datasets
- `res`: For datasets with multiple resolutions

Normal `GeoStack` keywords are passed to the constructor.
"""
GeoStack(T::Type{<:RasterDataSource}; kw...) = GeoStack(T, RDS.layers(T); kw...) 
GeoStack(T::Type{<:RasterDataSource}, layer::Symbol; kw...) = GeoStack(T, (layer,); kw...) 
function GeoStack(T::Type{<:RasterDataSource}, layers::LayerItr; kw...)
    rds_kw, gd_kw = _filterkw(kw)
    filenames = map(l -> getraster(T, l; rds_kw...), layers)
    GeoStack(filenames; keys=_layerkey(T, layers), gd_kw...)
end

"""
    GeoSeries(T::Type{<:RasterDataSource}, [layers::Union{Symbol,AbstractArray,Tuple}]; kw...) => AbstractGeoSeries

Load a `RasterDataSource` as an `AbstractGeoSeries`. `T`, `args` are
are passed to `getraster`, while `kw` args are for both `getraster` and
`AbstractGeoSeries`.

# Keywords

- `month`: For `Climate` datasets
- `date`: For `Weather` datasets
- `res`: For datasets with multiple resolutions

Normal `GeoStack` keywords are passed to the constructor.
"""
GeoSeries(T::Type{<:RasterDataSource}; kw...) = GeoSeries(T, RDS.layers(T); kw...) 
GeoSeries(T::Type{<:RasterDataSource}, layer::Symbol; kw...) = GeoSeries(T, (layer,); kw...) 
# Int month time-series
function GeoSeries(T::Type{WorldClim{Climate}}, layers::LayerItr;
    res=RDS.defres(T), month=1:12, 
    window=nothing, resize=nothing, crs=nothing, mappedcrs=nothing, kw...
)
    timedim = Ti(month; mode=Sampled(span=Regular(1), sampling=Intervals(Start())))
    stacks = [GeoStack(T, layers; res=res, month=m, window, crs, mappedcrs, resize) for m in month]
    GeoSeries(stacks, timedim; kw...)
end
# DateTime time-series
function GeoSeries(T::Type{<:Union{WorldClim{Weather},ALWB,AWAP}}, layers::LayerItr; 
    date, window=nothing, resize=nothing, crs=nothing, mappedcrs=nothing, kw...
)
    step = _seriesstep(T)
    dates = RDS._date_sequence(date, step)
    timedim = Ti(dates; mode=Sampled(Ordered(), Regular(step), Intervals(Start())))
    stacks = [stack(T, layers; date=d, window, crs, mappedcrs, resize) for d in dates]
    GeoSeries(stacks, timedim; kw...)
end

_sourcekw(T) = ()
_sourcekw(T::Type{AWAP}) = (crs=EPSG(4326),)

_layerkey(T::Type{<:RasterDataSource}, keys::LayerItr) = map(k -> _layerkey(T, k), keys) 
_layerkey(T::Type{<:Union{CHELSA{BioClim},WorldClim{BioClim}}}, key::Int) = Symbol(string("BIO", key))
_layerkey(T::Type{<:RasterDataSource}, key) = Symbol(key)

_seriesstep(T::Type{<:ALWB{M,P}}) where {M,P} = P(1)
_seriesstep(T::Type{<:WorldClim{<:Weather}}) = Month(1)
_seriesstep(T::Type{<:AWAP}) = Day(1)

function _filterkw(kw)
    rds = []; gd = []
    for p in kw
        dest = first(p) in (:date, :month, :res) ? rds : gd
        push!(dest, p)
    end
    rds, gd
end
