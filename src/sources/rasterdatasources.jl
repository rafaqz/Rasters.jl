using .RasterDataSources

using .RasterDataSources: RasterDataSource

const RDS = RasterDataSources

"""
    GeoArray(T::Type{<:RasterDataSource}, [layer]; kw...) => GeoArray

Load a `RasterDataSource` as an `GeoArray`. `T` and `layers` are are passed to
`RasterDataSources.getraster`, while `kw` args are for both `getraster` and `GeoArray`.

# Keywords

- `month`: an `Int` between `1` and `12`, usually for `Climate` datasets
- `date`: a `DateTime` object, usually for `Weather` datasets.
- `res`: a `String` resolution, for datasets with multiple resolutions.

Other `GeoArray` keywords are passed to the `GeoArray` constructor.

See the docs for 
[`RasterDatasources.getraster`](http://docs.ecojulia.org/RasterDataSources.jl/stable/#getraster)
for more specific details about data sources, layers and keyword arguments.
"""
function GeoArray(T::Type{<:RasterDataSource}, layer; crs=_source_crs(T), kw...)
    rds_kw, gd_kw = _filterkw(kw)
    filename = getraster(T, layer; rds_kw...)
    GeoArray(filename; name=RDS.layerkeys(T, layer), crs, gd_kw...)
end

"""
    GeoStack(T::Type{<:RasterDataSource}, [layers::Union{Symbol,AbstractArray,Tuple}]; kw...) => GeoStack

Load a `RasterDataSource` as an `GeoStack`. `T` and `layers` are passed to
`RasterDataSources.getraster`, while `kw` args are for both `getraster` and `GeoStack`.

# Keywords

- `month`: an `Int` between `1` and `12`, usually for `Climate` datasets.
- `date`: a `DateTime` object, usually for `Weather` datasets.
- `res`: a `String` resolution, for datasets with multiple resolutions.

Other `GeoStack` keywords are passed to the `GeoStack` constructor.

See the docs for 
[`RasterDatasources.getraster`](http://docs.ecojulia.org/RasterDataSources.jl/stable/#getraster)
for more specific details about data sources, layers and keyword arguments.
"""
GeoStack(T::Type{<:RasterDataSource}; kw...) = GeoStack(T, RDS.layers(T); kw...) 
GeoStack(T::Type{<:RasterDataSource}, layer::Symbol; kw...) = GeoStack(T, (layer,); kw...) 
function GeoStack(T::Type{<:RasterDataSource}, layers::Tuple; crs=_source_crs(T), kw...)
    rds_kw, gd_kw = _filterkw(kw)
    filenames = map(l -> getraster(T, l; rds_kw...), layers)
    GeoStack(filenames; keys=RDS.layerkeys(T, layers), crs, gd_kw...)
end

"""
    GeoSeries(T::Type{<:RasterDataSource}, [layers::Union{Symbol,AbstractArray,Tuple}]; kw...) => AbstractGeoSeries

Load a `RasterDataSource` as an `AbstractGeoSeries`. `T`, `args` are are passed to
`RasterDataSource.getraster`, while `kw` args are for both `getraster` and `GeoSeries`.

# Keywords

- `month`: a `Vector` or range of `Int` between `1` and `12`, usually for `Climate` datasets.
- `date`: a `Vector` of `DateTime` objects, usually for `Weather` datasets.
- `res`: a `String` resolution, for datasets with multiple resolutions.

Other `GeoSeries` keywords are passed to the `GeoSeries` constructor.

See the docs for 
[`RasterDatasources.getraster`](http://docs.ecojulia.org/RasterDataSources.jl/stable/#getraster)
for more specific details about data sources, layers and keyword arguments.
"""
GeoSeries(T::Type{<:RasterDataSource}; kw...) = GeoSeries(T, RDS.layers(T); kw...) 
# DateTime time-series
function GeoSeries(T::Type{<:RasterDataSource}, layers; 
    resize=_mayberesize(T), crs=_source_crs(T), mappedcrs=nothing, kw...
)
    monthdim = if haskey(values(kw), :month) values(kw)[:month] isa AbstractArray
        Dim{:month}(values(kw)[:month]; lookup=Sampled(; sampling=Intervals(Start())))
    else
        nothing
    end
    datedim = if haskey(values(kw), :date)
        dates = if values(kw)[:date] isa Tuple
            dates = RasterDataSources.date_sequence(T, values(kw)[:date])
            Ti(dates; lookup=Sampled(; sampling=Intervals(Start())))
        elseif values(kw)[:date] isa AbstractArray
            Ti(dates; lookup=Sampled(; sampling=Intervals(Start())))
        else
            nothing
        end
    else
        nothing
    end
    if isnothing(monthdim) && isnothing(datedim)
        throw(ArgumentError("A GeoSeries can only be constructed from a data source with `date` or `month` keywords that are AbstractArray. For other sources, use GeoStack or GeoArray directly"))
    end

    filenames = getraster(T, layers; kw...)
    can_duplicate = RDS.has_constant_dims(T) && RDS.has_constant_metadata(T)

    if filenames isa AbstractVector{<:AbstractVector}
        series = [GeoSeries(inner_fns, monthdim; resize, crs, mappedcrs, duplicate_first=can_duplicate) for inner_fns in filenames]
        return GeoSeries(series, datedim)
    elseif filenames isa AbstractVector
        dim = isnothing(datedim) ? monthdim : datedim
        return GeoSeries(filenames, dim; resize, crs, mappedcrs, duplicate_first=can_duplicate)
    end
end

_mayberesize(T) = RDS.has_matching_layer_size(T) ? nothing : crop

_source_crs(T) = nothing
_source_crs(T::Type{AWAP}) = crs=EPSG(4326)
_source_crs(T::Type{ALWB}) = crs=EPSG(4326)

function _filterkw(kw)
    rds = []; gd = []
    for p in kw
        dest = first(p) in (:date, :month, :res) ? rds : gd
        push!(dest, p)
    end
    rds, gd
end
