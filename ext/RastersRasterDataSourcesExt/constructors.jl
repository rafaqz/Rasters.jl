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
RA.Raster(T::Type{WorldClim{Elevation}}; kw...) = Raster(T, :elev; kw...)

# SRTM - special handling for tile-based data
"""
    Raster(::Type{SRTM}; bounds=nothing, tile_index=nothing, kw...) => Raster

Load SRTM elevation data as a `Raster`. Either `bounds` or `tile_index` must be specified.

# Keywords
- `bounds`: A tuple of (xmin, xmax, ymin, ymax) in WGS84 coordinates
- `tile_index`: A `CartesianIndex` or `CartesianIndices` specifying tile(s) to load
- `lazy`: If `true`, return a lazy raster (default `false`)

Other keywords are passed to the `Raster` constructor.

SRTM tiles are 5x5 degree blocks. When multiple tiles are requested, they are
lazily concatenated using `DiskArrays.ConcatDiskArray`. Missing ocean tiles
are filled with the `missingval`.
"""
function RA.Raster(T::Type{RDS.SRTM};
    bounds=nothing, tile_index=nothing, lazy=false, crs=EPSG(4326), kw...
)
    filenames = RDS.getraster(T; bounds, tile_index)

    raster = if filenames isa String
        # Single tile
        Raster(filenames; crs, lazy, kw...)
    else
        # Multiple tiles - load lazily and concatenate
        _load_srtm_tiles(filenames; crs, lazy, kw...)
    end
    # SRTM HGT files don't have embedded CRS metadata, so ensure CRS is set
    return RA.crs(raster) === nothing ? setcrs(raster, crs) : raster
end

function _load_srtm_tiles(filenames::AbstractMatrix; crs, lazy, kw...)
    # Find first valid tile to get reference dimensions/type
    first_valid_idx = findfirst(!ismissing, filenames)
    isnothing(first_valid_idx) && throw(ArgumentError("No valid SRTM tiles found for the specified region"))

    ref_raster = Raster(filenames[first_valid_idx]; crs, lazy=true, kw...)
    tile_size = size(ref_raster)
    T = eltype(ref_raster)
    mv = missingval(ref_raster)
    fill_val = ismissing(mv) || isnothing(mv) ? typemin(Int16) : mv

    # Build array of tile data, using Fill for missing ocean tiles
    tile_arrays = map(filenames) do f
        if ismissing(f)
            FillArrays.Fill(T(fill_val), tile_size)
        else
            parent(Raster(f; crs, lazy=true, kw...))
        end
    end

    # Concatenate: SRTM tiles are arranged with Y increasing downward in the matrix
    # Row 1 = northern tiles, last row = southern tiles
    # Columns go west to east
    # We need to concatenate columns along X, rows along Y
    nrows, ncols = size(tile_arrays)

    if ncols > 1
        # Concatenate each row along X (dimension 1 in the data)
        row_concats = [DA.ConcatDiskArray(tile_arrays[i, :]) for i in 1:nrows]
    else
        row_concats = [tile_arrays[i, 1] for i in 1:nrows]
    end

    if nrows > 1
        # Concatenate rows along Y (dimension 2 in the data)
        data = DA.ConcatDiskArray(row_concats)
    else
        data = row_concats[1]
    end

    # Build combined dimensions
    combined_dims = _srtm_combined_dims(ref_raster, filenames)

    result = rebuild(ref_raster; data, dims=combined_dims)
    return lazy ? result : read(result)
end

function _srtm_combined_dims(ref_raster, filenames::AbstractMatrix)
    # Each tile spans 5 degrees
    # The filenames matrix has rows for Y tiles (north to south) and cols for X tiles (west to east)
    nrows, ncols = size(filenames)

    ref_x = dims(ref_raster, X)
    ref_y = dims(ref_raster, Y)
    tile_nx = length(ref_x)
    tile_ny = length(ref_y)

    # Calculate the step size from the reference tile
    x_step = step(ref_x)
    y_step = step(ref_y)

    # Build combined X dimension (concatenate columns)
    # Rebuild the lookup with new data to preserve Lookup structure (including CRS)
    x_start = first(ref_x)
    x_len = tile_nx * ncols
    new_x_data = LinRange(x_start, x_start + x_step * (x_len - 1), x_len)
    combined_x = rebuild(ref_x; val=rebuild(lookup(ref_x); data=new_x_data))

    # Build combined Y dimension (concatenate rows)
    # Y goes from north (first row) to south (last row) in SRTM
    y_start = first(ref_y)
    y_len = tile_ny * nrows
    new_y_data = LinRange(y_start, y_start + y_step * (y_len - 1), y_len)
    combined_y = rebuild(ref_y; val=rebuild(lookup(ref_y); data=new_y_data))

    return (combined_x, combined_y)
end
# WorldClim Future BioClim - all layers in one file, need to select by band
function RA.Raster(T::Type{<:WorldClim{<:Future{BioClim}}};
    crs=_source_crs(T), kw...
)
    rds_kw, ra_kw = _filterkw(T, kw)
    filename = getraster(T; rds_kw...)
    Raster(filename; crs, ra_kw...)
end
function RA.Raster(T::Type{<:WorldClim{<:Future{BioClim}}}, layer::Union{Symbol,Int};
    crs=_source_crs(T), lazy=false, kw...
)
    rds_kw, ra_kw = _filterkw(T, kw)
    filename = getraster(T, layer; rds_kw...)
    ras_all = Raster(filename; name=RDS.bioclim_key(layer), crs, ra_kw...)
    ras = view(ras_all, Band(RDS.bioclim_int(layer)))
    return lazy ? ras : read(ras)
end

# MODIS - requires specific keywords
function RA.Raster(T::Type{<:RDS.MODIS}, layer::Union{Symbol,Int}=first(RDS.layerkeys(T));
    lat, lon, km_ab, km_lr, date, crs=_source_crs(T), kw...
)
    rds_kw = (; lat, lon, km_ab, km_lr, date)
    filename = RDS.getraster(T, layer; rds_kw...)
    if filename isa AbstractVector
        # Multiple dates - return series
        dates = date isa Tuple ? RDS.date_sequence(T, date; lat, lon) : date
        return RasterSeries(filename, Ti(dates); crs, kw...)
    end
    Raster(filename; name=RDS.layerkeys(T, layer), crs, kw...)
end

function RA.RasterStack(T::Type{<:RDS.MODIS}, layers::Tuple=RDS.layerkeys(T);
    lat, lon, km_ab, km_lr, date, crs=_source_crs(T), kw...
)
    rds_kw = (; lat, lon, km_ab, km_lr, date)
    filenames = map(l -> RDS.getraster(T, l; rds_kw...), layers)
    RasterStack(filenames; name=RDS.layerkeys(T, layers), crs, kw...)
end

function RA.RasterStack(T::Type{<:RDS.MODIS}, layer::Symbol;
    lat, lon, km_ab, km_lr, date, kw...
)
    RasterStack(T, (layer,); lat, lon, km_ab, km_lr, date, kw...)
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
# WorldClim Future BioClim - all layers in one .tif file with Band dimension
function RA.RasterStack(T::Type{<:WorldClim{<:Future{BioClim}}}, layers::Tuple; crs=_source_crs(T), lazy=false, kw...)
    rds_kw, ra_kw = _filterkw(T, kw)
    keys = map(RDS.bioclim_key, layers)
    filename = RDS.getraster(T, layers; rds_kw...)
    raster = Raster(filename; lazy=true, ra_kw...)
    stack = RasterStack(eachslice(raster; dims=Band)...; name=RDS.layerkeys(T), crs)[keys]
    return lazy ? stack : read(stack)
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
