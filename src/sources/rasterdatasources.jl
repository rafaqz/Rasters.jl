using .RasterDataSources

using .RasterDataSources: RasterDataSource

const RDS = RasterDataSources

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
function GeoArray(T::Type{<:RasterDataSource}, layer; crs=_source_crs(T), kw...)
    rds_kw, gd_kw = _filterkw(kw)
    filename = getraster(T, layer; rds_kw...)
    GeoArray(filename; name=RDS.layerkeys(T, layer), crs, gd_kw...)
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
function GeoStack(T::Type{<:RasterDataSource}, layers::Tuple; crs=_source_crs(T), kw...)
    rds_kw, gd_kw = _filterkw(kw)
    filenames = map(l -> getraster(T, l; rds_kw...), layers)
    GeoStack(filenames; keys=RDS.layerkeys(T, layers), crs, gd_kw...)
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
# DateTime time-series
function GeoSeries(T::Type{<:RasterDataSource}, layers; 
    resize=_mayberesize(T), crs=_source_crs(T), mappedcrs=nothing, kw...
)
    monthdim = if haskey(kw.data, :month) kw.data[:month] isa AbstractArray
        Dim{:month}(kw.data[:month]; mode=Sampled(; sampling=Intervals(Start())))
    else
        nothing
    end
    datedim = if haskey(kw.data, :date)
        dates = if kw.data[:date] isa Tuple
            dates = RasterDataSources.date_sequence(T, kw.data[:date])
            Ti(dates; mode=Sampled(; sampling=Intervals(Start())))
        elseif kw.data[:date] isa AbstractArray
            Ti(dates; mode=Sampled(; sampling=Intervals(Start())))
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
