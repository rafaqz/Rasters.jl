struct _Undefined end
struct _Defined end

"""
    rasterize(obj; to, fill, kw...)

Rasterize a GeoInterface.jl compatable geometry or feature,
or a Tables.jl table with a `:geometry` column of GeoInterface.jl objects,
or `X`, `Y` points columns.

# Arguments

- `obj`: a GeoInterface.jl `AbstractGeometry`, or a nested `Vector` of `AbstractGeometry`,
    or a Tables.jl compatible object containing points and values columns or an `:geometry`.

# Keywords

These are detected automatically from `obj` where possible.

- `fill`: the value to fill a polygon with. A `Symbol` or tuple of `Symbol` will
    be used to retrieve properties from features or column values from table rows.
- `atol`: an absolute tolerance for rasterizing to dimensions with `Points` sampling.
- `filename`: a filename to write to directly, useful for large files.
- `suffix`: a string or value to append to the filename.
    A tuple of `suffix` will be applied to stack layers. `keys(st)` are the default.

## Geometry keywords

These can be used when a `GeoInterface.AbstractGeometry` is passed in.

$GEOM_KEYWORDS

# Example

Rasterize a shapefile for China and plot, with a border.

```jldoctest
using Rasters, Plots, Dates, Shapefile, Downloads
using Rasters.LookupArrays

# Download a borders shapefile
shapefile_url = "https://github.com/nvkelso/natural-earth-vector/raw/master/10m_cultural/ne_10m_admin_0_countries.shp"
shapefile_name = "country_borders.shp"
isfile(shapefile_name) || Downloads.download(shapefile_url, shapefile_name)

# Loade the shapes for china
china_border = Shapefile.Handle(shapefile_name).shapes[10]

# Rasterize the border polygon
china = rasterize(china_border; res=0.1, missingval=0, fill=1, boundary=:touches)

# And plot
p = plot(china; color=:spring)
plot!(p, china_border; fillalpha=0, linewidth=0.6)

savefig("build/china_rasterized.png"); nothing

# output

```

![rasterize](china_rasterized.png)

$EXPERIMENTAL
"""
function rasterize(data; to=nothing, fill, kw...)
    return _rasterize(to, data; fill, kw...)
end
function _rasterize(to::AbstractRaster, data;
    missingval=missingval(to), name=name(to), kw...
)
    _rasterize(dims(to), data; missingval, name, kw...)
end
function _rasterize(to::AbstractRasterStack, data; fill, name=keys(to), kw...)
    _rasterize(dims(to), data; fill, name=_filter_name(name, fill), kw...)
end
function _rasterize(to::Nothing, data; fill, kw...)
    to = _extent(data)
    _rasterize(to, data; fill, name, kw...)
end
function _rasterize(to::Extents.Extent{K}, data;
    fill, name=_filter_name(nothing, fill),
    res::Union{Nothing,Real,NTuple{<:Any,<:Real}}=nothing,
    size::Union{Nothing,Int,NTuple{<:Any,Int}}=nothing,
    kw...
) where K
    to_dims = _extent2dims(to; size, res, kw...)
    return _rasterize(to_dims, data; fill, name, kw...)
end
function _rasterize(to::DimTuple, data; fill, name=_filter_name(nothing, fill), kw...)
    _rasterize(to, GeoInterface.trait(data), data; fill, name, kw...)
end
function _rasterize(to::DimTuple, ::GI.AbstractFeatureCollectionTrait, fc; name, fill, kw...)
    # TODO: how to handle when there are fillvals with different types
    fillval = _featurefillval(GI.getfeature(fc, 1), fill)
    name = _filter_name(name, fill)
    return _create_rasterize_dest(fillval, to; name, kw...) do dest
        rasterize!(dest, fc; fill, kw...)
    end
end
function _rasterize(to::DimTuple, ::GI.AbstractFeatureTrait, feature; fill, name, kw...)
    fillval = _featurefillval(feature, fill)
    name = _filter_name(name, fill)
    return _create_rasterize_dest(fillval, to; name, kw...) do dest
        rasterize!(dest, feature; fill, kw...)
    end
end
function _rasterize(to::DimTuple, ::GI.AbstractGeometryTrait, geom; fill, kw...)
    return _create_rasterize_dest(fill, to; kw...) do dest
        rasterize!(dest, geom; fill, kw...)
    end
end
function _rasterize(to::DimTuple, ::Nothing, data; fill, name, kw...)
    name = _filter_name(name, fill)
    if Tables.istable(data)
        schema = Tables.schema(data)
        fillval = if isnothing(fill)
            throw(ArgumentError("`fill` must be a value or table column name or names"))
        elseif fill isa Symbol
            zero(Tables.columntype(schema, fill))
        elseif fill isa NTuple{<:Any,Symbol}
            map(n -> zero(Tables.columntype(schema, n)), fill)
        else
            fill
        end
        return _create_rasterize_dest(fillval, to; name, kw...) do dest
            rasterize!(dest, data; fill, kw...)
        end
    else
        return _create_rasterize_dest(fill, to; name, kw...) do dest
            rasterize!(dest, data; fill, kw...)
        end
    end
end

"""
    rasterize!(dest, data; fill, atol)

Rasterize the geometries in `data` into the [`Raster`](@ref) or [`RasterStack`](@ref) `dest`,
using the values specified by `fill`.

# Arguments

- `dest`: a `Raster` or `RasterStack` to rasterize into.
- `data`: a GeoInterface.jl compatible object or an `AbstractVector` of such objects,
    or a Tables.jl compatible table containing `:geometry` column of GeoInterface compatible objects or
    columns with point names `X` and `Y`.
- `fill`: the value to fill a polygon with. A `Symbol` or tuple of `Symbol` will
    be used to retrieve properties from features or column values from table rows.

# Keywords

These are detected automatically from `A` and `data` where possible.

- `atol`: an absolute tolerance for rasterizing points to dimensions with `Points` sampling.

## Geometry keywords

These can be used when a `GeoInterface.AbstractGeometry` is passed in.

- `shape`: Force `data` to be treated as `:polygon`, `:line` or `:point`.

And specifically for `shape=:polygon`:

- `boundary`: include pixels where the `:center` is inside the polygon, where
    the polygon `:touches` the pixel, or that are completely `:inside` the polygon.

## Table keywords

- `name`: A `Symbol` to return a `Raster` from a single column,
    or `Tuple` of `Symbol` to return a `RasterStack` from multiple columns.

# Example

```jldoctest
using Rasters, Plots, Dates, Shapefile, GeoInterface, Downloads
using Rasters.LookupArrays

# Download a borders shapefile
shapefile_url = "https://github.com/nvkelso/natural-earth-vector/raw/master/10m_cultural/ne_10m_admin_0_countries.shp"
shapefile_name = "country_borders.shp"
isfile(shapefile_name) || Downloads.download(shapefile_url, shapefile_name)

# Load the shapes for indonesia
indonesia_border = Shapefile.Handle(shapefile_name).shapes[1]

# Make an empty EPSG 4326 projected Raster of the area of Indonesia
dimz = X(90.0:0.1:145; mode=Projected(; sampling=Points(), crs=EPSG(4326))),
       Y(-15.0:0.1:10.9; mode=Projected(; sampling=Points(), crs=EPSG(4326)))

A = Raster(zeros(UInt32, dimz); missingval=UInt32(0))

# Rasterize each island with a different number
for (i, shp) in enumerate(GeoInterface.getring(indonesia_border))
    rasterize!(A, shp; fill=i)
end

# And plot
p = plot(Rasters.trim(A); color=:spring)
plot!(p, indonesia_border; fillalpha=0, linewidth=0.7)

savefig("build/indonesia_rasterized.png"); nothing

# output

```

![rasterize](indonesia_rasterized.png)

$EXPERIMENTAL
"""
rasterize!(x::RasterStackOrArray, data; fill, kw...) =
    _rasterize!(x, GI.trait(data), data; fill, kw...)

function _rasterize!(x, ::GI.AbstractFeatureCollectionTrait, fc; fill, kw...)
    function _rasterize_feature_inner(x, fc, fillkey; kw...)
        for feature in GI.getfeature(fc)
            geom = GI.geometry(feature)
            rasterize!(x, geom; fill=_featurefillval(feature, fillkey), kw...)
        end
    end
    buffer = _fillraster(commondims(x, (XDim, YDim)), Bool; missingval=false)
    if fill isa Union{Symbol,NTuple{<:Any,Symbol}}
        # Rasterize features separately: the fill may change per feature
        # Lift key Symbol to a type to avoid runtime lookups for every point.
        fillkey = Val{fill}()
        # Use a function barrier
        _rasterize_feature_inner(x, fc, fillkey; _buffer=buffer, kw...)
    else
        # Rasterize all features together: the fill is the same
        boolmask!(buffer, fc; kw...)
        _fill!(x, buffer, fill, _Defined())
    end
    return x
end
function _rasterize!(x, ::GI.AbstractFeatureTrait, feature; fill, kw...)
    rasterize!(x, GI.geometry(feature); fill=_featurefillval(feature, fill), kw...)
end
function _rasterize!(x, ::GI.AbstractGeometryTrait, geom; fill, _buffer=nothing, kw...)
    # TODO fix DimensionalData selectors so this works without _pad
    ranges = _pad(size(x), Base.to_indices(x, DD.dims2indices(x, _extent(geom))))
    x1 = view(x, ranges...)
    length(x1) > 0 || return x
    _buffer = if isnothing(_buffer)
        _fillraster(commondims(x1, (XDim, YDim)), Bool; missingval=false)
    else
        view(_buffer, ranges...)
    end
    boolmask!(_buffer, geom; kw...)
    _fill!(x1, _buffer, fill, _Defined())
    return x
end
function _rasterize!(x, trait::GI.AbstractPointTrait, point; fill, kw...)
    _fill_point!(x, trait, point; fill, kw...)
    return x
end
function _rasterize!(x, trait::Nothing, data; fill, kw...)
    _buffer = _fillraster(commondims(x, (XDim, YDim)), Bool; missingval=false)

    if Tables.istable(data)
        schema = Tables.schema(data)
        geomcolname = first(GI.geometrycolumns(data))
        cols = Tables.columns(data)
        fillcol = _fillorcol(data, fill)
        if geomcolname in Tables.columnnames(cols)
            geomcol = Tables.getcolumn(data, geomcolname)
            _geom_table_inner(x, geomcol, fillcol; _buffer, kw...)
        else
            pointkeys = reduce(DEFAULT_TABLE_DIM_KEYS; init=()) do acc, key
                key in schema.names ? (acc..., key) : acc
            end
            pointcols = map(k -> Tables.getcolumn(data, k), pointkeys)
            _point_table_inner(x, pointcols, fillcol; _buffer, kw...)
        end
    elseif data isa AbstractArray
        for geom in data
            rasterize!(x, geom; fill, _buffer, kw...)
        end
    else
        throw(ArgumentError("data should be either a GeoInterface.jl compatible object, or a table with a :geometry column"))
    end
    return x
end

function _geom_table_inner(x, geomcol, fill; kw...)
    for i in eachindex(geomcol)
        fillval = _fillval(fill, i)
        geom = geomcol[i]
        rasterize!(x, geom; fill=fillval, kw...)
    end
end

function _point_table_inner(x, pointcols, fill; kw...)
    for i in eachindex(first(pointcols))
        fillval = _fillval(fill, i)
        point = map(col -> col[i], pointcols)
        _fill_point!(x, point; fill=fillval, kw...)
    end
end

# Utils

_fillorcol(data, fill::NTuple{<:Any,Symbol}) = map(f -> _fillorcol(data, f), fill)
_fillorcol(data, fill::Symbol) = Tables.getcolumn(data, fill)
_fillorcol(data, fill) = fill

# _fillval
# Get fill value from a table row, or use fill itself
_fillval(fill::Tuple, i) = map(f -> _fillval(f, i), fill)
_fillval(fillcol::AbstractArray, i) = fillcol[i]
_fillval(fill, i) = fill

# _featurefillval
# Get fill value from a feature, or use fill itself
_featurefillval(feature, fill::Nothing) = first(GI.properties(feature))
_featurefillval(feature, fill::Symbol) = GI.properties(feature)[fill]
_featurefillval(feature, fill::Val) = _featurefillval(feature, _unwrap(fill))
function _featurefillval(feature, fill::NTuple{<:Any,Symbol})
    map(fill) do key
        getproperty(GI.properties(feature), key)
    end
end
_featurefillval(feature, fill) = fill

function _pad(size, ranges)
    map(ranges, size) do range, s
        max(1, first(range) - 1):min(s, last(range) + 1)
    end
end

function _create_rasterize_dest(f, fill, dims; name=nothing, kw...)
    _create_rasterize_dest(f, fill, name, dims; kw...)
end
function _create_rasterize_dest(f, fill::Union{Tuple,NamedTuple}, keys, dims; kw...)
    _create_rasterize_dest(f, fill, DD.uniquekeys(fill), dims; kw...)
end
function _create_rasterize_dest(f, fill::Union{Tuple,NamedTuple}, keys::Union{Tuple,NamedTuple}, dims;
    filename=nothing, missingval=nothing, metadata=NoMetadata(), suffix=nothing, kw...
)
    layers = map(keys, values(fill)) do key, val
        missingval = isnothing(missingval) ? _writeable_missing(filename, typeof(val)) : missingval
        _alloc_rasterize(filename, val, dims; name, metadata, missingval, suffix=key) do a
            a .= missingval
        end
    end
    st = RasterStack(layers, dims; keys, metadata)
    open(f, RasterStack(layers, dims; keys, metadata))
    return st
end
function _create_rasterize_dest(f, fill, name, dims;
    filename=nothing, missingval=nothing, metadata=NoMetadata(), suffix=nothing, kw...
)
    missingval = isnothing(missingval) ? _writeable_missing(filename, typeof(val)) : missingval
    A = _alloc_rasterize(filename, fill, dims; name, metadata, missingval, suffix) do a
        a .= missingval
        f(a)
    end
    return Raster(A, dims; name, missingval, metadata)
end

function _alloc_rasterize(f, filename, fill, to; missingval, kw...)
    T = fill isa AbstractArray ? eltype(fill) : typeof(fill)
    T1 = promote_type(typeof(missingval), T)
    A = create(filename, T1, to; missingval, kw...)
    # TODO f should apply to the file when it is initially created
    # instead of reopening but we need a `create(f, filename, ...)` method
    open(A; write=true) do A
        A .= missingval
        f(A)
    end
    return A
end

function _fill!(st::AbstractRasterStack, B, fill, args...)
    map((a, f) -> _fill!(a, B, f, args...), values(st), fill)
    return st
end
# If the array is initialised, we can use the existing values
function _fill!(A::AbstractRaster{T}, B, fill, init::_Defined, missingval=nothing) where T
    broadcast_dims!(A, A, B) do a, b
        val = b ? (fill isa Function ? fill(a) : fill) : a
        convert(T, val) # In case we are writing to disk
    end
end
# If the array is not yet initialised, we have to fill with fill and misssingval
function _fill!(A::AbstractRaster{T}, B, fill, init::_Undefined, missingval) where T
    fill = convert(T, fill) # In case we are writing to disk
    missingval = convert(T, missingval)
    broadcast_dims!(A, B) do b
        b ? fill : missingval
    end
end

function _at_or_contains(d, v, atol)
    selector = sampling(d) isa Intervals ? Contains(v) : At(v; atol=atol)
    DD.basetypeof(d)(selector)
end

_filter_name(name, fill::NamedTuple) = keys(fill)
_filter_name(name::NamedTuple, fill::NamedTuple) = keys(fill)
_filter_name(name::Nothing, fill::Nothing) = nothing
_filter_name(name::DimensionalData.NoName, fill::Nothing) = nothing
_filter_name(name::Union{NamedTuple,Tuple,Array}, fill::NTuple{<:Any,Symbol}) = fill
function _filter_name(name::Union{NamedTuple,Tuple,Array}, fill::Union{Tuple,Array})
    length(name) == length(fill) || throw(ArgumentError("`name` keyword (possibly from `to` object) does not match length of fill. A fix is to use a `NamedTuple` for `fill`."))
    return name isa NamedTuple ? keys(name) : name
end
function _filter_name(name, fill)
    fill isa Union{Symbol,NTuple{<:Any,Symbol}} ? fill : name
end
