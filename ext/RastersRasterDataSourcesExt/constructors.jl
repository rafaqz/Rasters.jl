"""
    Raster(T::Type{<:RasterDataSource}, [layer]; kw...) => Raster

Load a `RasterDataSource` as an `Raster`. `T` and `layers` are are passed to
`RasterDataSources.getraster`, while `kw` args are for both `getraster` and `Raster`.

# Keywords

- `month`: an `Int` between `1` and `12`, usually for `Climate` datasets
- `date`: a `DateTime` object, usually for `Weather` datasets.
- `res`: a `String` resolution, for datasets with multiple resolutions.

`MODIS` datasets require a specific set of keyword arguments:

- `lat` and `lon`: Coordinates in decimal degrees (`Float`s ) of the center of the raster
- `km_ab` and `km_lr`: Kilometers above/below and left/right (`Integer`s up to 100) of the center of the raster

Other `Raster` keywords are passed to the `Raster` constructor.

See the docs for 
[`RasterDatasources.getraster`](http://docs.ecojulia.org/RasterDataSources.jl/stable/#getraster)
for more specific details about data sources, layers and keyword arguments.
"""
function RA.Raster(T::Type{<:RDS.RasterDataSource}, layer::Union{Symbol,Int}; 
    crs=_source_crs(T), kw...
)
    rds_kw, ra_kw = _filterkw(T, kw)
    filename = getraster(T, layer; rds_kw...)
    Raster(filename; name=RDS.layerkeys(T, layer), crs, ra_kw...)
end
function RA.Raster(T::Type{<:WorldClim{<:Future{BioClim, CMIP6}}}; 
    crs=_source_crs(T), kw...
)
    rds_kw, ra_kw = _filterkw(T, kw)
    filename = getraster(T; rds_kw...)
    Raster(filename; crs, ra_kw...)
end
function RA.Raster(T::Type{<:WorldClim{<:Future{BioClim, CMIP6}}}, layer::Union{Symbol,Int}; 
    crs=_source_crs(T), lazy = false, kw...
)
    rds_kw, ra_kw = _filterkw(T, kw)
    filename = getraster(T, layer; rds_kw...)
    ras_all = Raster(filename; name=RDS.bioclim_key(layer), crs, ra_kw...)
    ras = view(ras_all, Band(RDS.bioclim_int(layer)))
    if ~lazy
        read(ras)
    end
    return ras
end
"""
    RasterStack(T::Type{<:RasterDataSource}, [layers::Union{Symbol,AbstractArray,Tuple}]; kw...) => RasterStack

Load a `RasterDataSource` as an `RasterStack`. `T` and `layers` are passed to
`RasterDataSources.getraster`, while `kw` args are for both `getraster` and `RasterStack`.

# Keywords

- `month`: an `Int` between `1` and `12`, usually for `Climate` datasets.
- `date`: a `DateTime` object, usually for `Weather` datasets.
- `res`: a `String` resolution, for datasets with multiple resolutions.

`MODIS` datasets require a specific set of keyword arguments:

- `lat` and `lon`: Coordinates in decimal degrees (`Float`s ) of the center of the raster
- `km_ab` and `km_lr`: Kilometers above/below and left/right (`Integer`s up to 100) of the center of the raster

Other `RasterStack` keywords are passed to the `RasterStack` constructor.

See the docs for 
[`RasterDatasources.getraster`](http://docs.ecojulia.org/RasterDataSources.jl/stable/#getraster)
for more specific details about data sources, layers and keyword arguments.
"""
RA.RasterStack(T::Type{<:RDS.RasterDataSource}; kw...) = RasterStack(T, RDS.layers(T); kw...) 
RA.RasterStack(T::Type{<:RDS.RasterDataSource}, layer::Symbol; kw...) = RasterStack(T, (layer,); kw...) 
function RA.RasterStack(T::Type{<:RDS.RasterDataSource}, layers::Tuple; crs=_source_crs(T), kw...)
    rds_kw, ra_kw = _filterkw(T, kw)
    filenames = map(l -> RDS.getraster(T, l; rds_kw...), layers)
    RasterStack(filenames; name=RDS.layerkeys(T, layers), crs, ra_kw...)
end
# WorldClim future has all layers in one .tif file - format this as a rasterdatastack in a dedicated dispatch
function RA.RasterStack(T::Type{<:WorldClim{<:Future{BioClim, CMIP6}}}, layers::Tuple; crs=_source_crs(T), lazy = false, kw...)
    rds_kw, ra_kw = _filterkw(T, kw)
    keys = map(RDS.bioclim_key, layers)
    filename = RDS.getraster(T, layers; rds_kw...)
    raster = Raster(filename; lazy = true, ra_kw...)
    stack = RasterStack(eachslice(raster; dims = Band)...; name=RDS.layerkeys(T), crs)[keys]
    if ~lazy
        stack = read(stack)
    end
    return stack
end

"""
    RasterSeries(T::Type{<:RasterDataSource}, [layers::Union{Symbol,AbstractArray,Tuple}]; kw...) => AbstractRasterSeries

Load a `RasterDataSource` as an `AbstractRasterSeries`. `T`, `args` are are passed to
`RasterDataSource.getraster`, while `kw` args are for both `getraster` and `RasterSeries`.

# Keywords

- `month`: a `Vector` or range of `Int` between `1` and `12`, usually for `Climate` datasets.
- `date`: a `Vector` of `DateTime` objects, usually for `Weather` datasets.
- `res`: a `String` resolution, for datasets with multiple resolutions.

`MODIS` datasets require a specific set of keyword arguments:

- `lat` and `lon`: Coordinates in decimal degrees (`Float`s ) of the center of the raster
- `km_ab` and `km_lr`: Kilometers above/below and left/right (`Integer`s up to 100) of the center of the raster

Other `RasterSeries` keywords are passed to the `RasterSeries` constructor.

See the docs for 
[`RasterDatasources.getraster`](http://docs.ecojulia.org/RasterDataSources.jl/stable/#getraster)
for more specific details about data sources, layers and keyword arguments.
"""
RA.RasterSeries(T::Type{<:RDS.RasterDataSource}; kw...) = RasterSeries(T, RDS.layers(T); kw...) 
# DateTime time-series
function RA.RasterSeries(T::Type{<:RDS.RasterDataSource}, layers; 
    resize=_mayberesize(T), 
    crs=_source_crs(T), 
    mappedcrs=RA.nokw, 
    kw...
)
    monthdim = if haskey(values(kw), :month) values(kw)[:month] isa AbstractArray
        Dim{:month}(values(kw)[:month]; lookup=Sampled(; sampling=Intervals(Start())))
    else
        nothing
    end
    datedim = if haskey(values(kw), :date)
        if values(kw)[:date] isa Tuple
            dates = RasterDataSources.date_sequence(T, values(kw)[:date]; kw...)
            Ti(dates; lookup=Sampled(; sampling=Intervals(Start())))
        elseif values(kw)[:date] isa AbstractArray
            dates = values(kw)[:date]
            Ti(dates; lookup=Sampled(; sampling=Intervals(Start())))
        else
            nothing
        end
    else
        nothing
    end
    if isnothing(monthdim) && isnothing(datedim)
        throw(ArgumentError("A RasterSeries can only be constructed from a data source with `date` or `month` keywords that are AbstractArray or Tuple. For other sources, use RasterStack or Raster directly"))
    end

    rds_kw, ra_kw = _filterkw(T, kw)
    filenames = RDS.getraster(T, layers; rds_kw...)
    can_duplicate = RDS.has_constant_dims(T) && RDS.has_constant_metadata(T)

    if filenames isa AbstractVector{<:AbstractVector}
        series = [RasterSeries(inner_fns, monthdim; resize, crs, mappedcrs, duplicate_first=can_duplicate, ra_kw...) for inner_fns in filenames]
        return RasterSeries(series, datedim)
    elseif filenames isa AbstractVector
        dim = isnothing(datedim) ? monthdim : datedim
        return RasterSeries(filenames, dim; resize, crs, mappedcrs, duplicate_first=can_duplicate, ra_kw...)
    else
        error("Returned filenames for `RasterSeries` must be a `Vector`. Got $filenames")
    end
end

_mayberesize(T) = RDS.has_matching_layer_size(T) ? RA.nokw : crop

_source_crs(T) = RA.nokw
_source_crs(T::Type{AWAP}) = crs=EPSG(4326)
_source_crs(T::Type{ALWB}) = crs=EPSG(4326)

function _filterkw(T, kw)
    rds = []; gd = []
    for p in kw
        dest = first(p) in RDS.getraster_keywords(T) ? rds : gd
        push!(dest, p)
    end
    rds, gd
end
