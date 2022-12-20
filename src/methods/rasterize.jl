"""
    rasterize(obj; to, fill, kw...)

Rasterize a GeoInterface.jl compatable geometry or feature,
or a Tables.jl table with a `:geometry` column of GeoInterface.jl objects,
or `X`, `Y` points columns.

# Arguments

- `obj`: a GeoInterface.jl `AbstractGeometry`, or a nested `Vector` of `AbstractGeometry`,
    or a Tables.jl compatible object containing a `:geometry` column or points and values columns.

# Keywords

These are detected automatically from `obj` where possible.

- `fill`: the value to fill a polygon with. A `Symbol` or tuple of `Symbol` will
    be used to retrieve properties from features or column values from table rows.
- `reduce`: a function to reduce the values of all geometries that cover or touch a
    pixel down to a single value. Default is `last`.
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
function rasterize end
rasterize(reduce::Function, data; kw...) = rasterize(data; reduce, kw...)
function rasterize(data; to=nothing, fill, kw...)
    return _rasterize(to, data; fill, kw...)
end

function _rasterize(to::AbstractRaster, data;
    missingval=missingval(to), name=name(to), kw...
)
    _rasterize(dims(to), data; missingval, name, kw...)
end
function _rasterize(to::AbstractRasterStack, data; fill, name=keys(to), kw...)
    _rasterize(dims(to), data; fill, name, kw...)
end
function _rasterize(to::Nothing, data; fill, kw...)
    to = _extent(data)
    _rasterize(to, data; fill, name, kw...)
end
function _rasterize(to::Extents.Extent{K}, data;
    fill, 
    res::Union{Nothing,Real,NTuple{<:Any,<:Real}}=nothing,
    size::Union{Nothing,Int,NTuple{<:Any,Int}}=nothing,
    kw...
) where K
    to_dims = _extent2dims(to; size, res, kw...)
    return _rasterize(to_dims, data; fill, kw...)
end
function _rasterize(to::DimTuple, data; fill, name=nothing, kw...)
    if Tables.istable(data)
        schema = Tables.schema(data)
        colnames = Tables.columnnames(Tables.columns(data))
        # TODO integrate this with _iterable_fill
        fillval = if fill isa Symbol
            fill in colnames || _fill_key_error(colnames, f)
            zero(Tables.columntype(schema, fill))
        elseif fill isa Tuple{Symbol,Vararg}
            map(fill) do f
                f in colnames || _fill_key_error(colnames, f)
                zero(Tables.columntype(schema, f)) 
            end
        else
            fill
        end
        name = _filter_name(name, fill)
        return _create_rasterize_dest(fillval, to; name, kw...) do dest
            rasterize!(dest, data; fill, kw...)
        end
    else
        _rasterize(to, GeoInterface.trait(data), data; fill, name, kw...)
    end
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
function _rasterize(to::DimTuple, ::Nothing, data; fill, kw...)
    # treat data as some kind of interable of geometries
    return _create_rasterize_dest(fill, to; kw...) do dest
        rasterize!(dest, data; fill, kw...)
    end
end


_fill_key_error(names, fill) = throw(ArgumentError("fill key $fill not found in table, use one of: $(Tuple(names))"))

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
function rasterize! end
rasterize!(reduce::Function, x::RasterStackOrArray, data; kw...) =
    rasterize!(x::RasterStackOrArray, data; reduce, kw...)
function rasterize!(x::RasterStackOrArray, data; fill, reduce=last, kw...)
    function _iterable_fill(data, fill::Symbol) 
        names = Tables.columnnames(Tables.columns(data))
        fill in names || _fill_key_error(names, fill)
        Tables.getcolumn(data, fill)
    end
    _iterable_fill(data, fill::Tuple{Symbol,Vararg}) =
        map(f -> _iterable_fill(data, f), fill)
    _iterable_fill(data, fill) = Iterators.cycle(fill)

    # Check if `data` is a Tables.jl compatible object
    if Tables.istable(data)
        schema = Tables.schema(data)
        geomcolname = first(GI.geometrycolumns(data))
        cols = Tables.columns(data)
        fill_itr = _iterable_fill(cols, fill)
        if geomcolname in Tables.columnnames(cols)
            geomcol = Tables.getcolumn(cols, geomcolname)
            _reduce_fill!(reduce, x, geoms, fill_itr; kw...)
        else
            dimscols = _auto_dim_columns(data, dims(x))
            pointcols = map(k -> Tables.getcolumn(data, k), map(DD.dim2key, dimscols))
            reduce == last || throw(ArgumentError("Can only reduce with `last` on point tables. Make a github issue at the Rasters.jl repository if you need this."))
            _buffer = _init_bools(commondims(x, (XDim, YDim)), Bool; missingval=false)
            _rasterize_point_table_inner!(x, pointcols, fill_itr; _buffer, kw...)
        end
    else
        # Otherwise treat as a GeoInterface compatible geometry
        _rasterize!(x, GI.trait(data), data; fill, kw...)
    end
end

function _rasterize_point_table_inner!(x, pointcols, fill_itr; kw...)
    for (point, fill) in  zip(zip(pointcols...), fill_itr)
        _fill_point!(x, point; fill, kw...)
    end
end

function _rasterize!(x, ::GI.AbstractFeatureCollectionTrait, fc; fill, kw...)
    function _iterable_fill(fc, keyorfill::Val)
        if _unwrap(keyorfill) isa Tuple
            (NamedTuple{_unwrap(key)}(map(p -> getproperty(GI.properties(f), p), _unwrap(key))) for f in GI.getfeature(fc))
        else
            (getproperty(GI.properties(f), _unwrap(key)) for f in GI.getfeature(fc))
        end
    end
    _iterable_fill(fc, keyorfill) = Iterators.cycle(keyorfill)

    if fill isa Union{Symbol,Tuple{Symbol,Vararg}}
        # Rasterize features separately: the fill may change per feature
        # Lift key Symbol to a type to avoid runtime lookups for every point.
        keyorfill = _iscolumnfill(fill) ? Val{fill}() : fill
        fill_itr = _iterable_fill(fc, keyorfill)
        geoms = (GI.geometry(feature) for feature in GI.getfeature(fc))
        _reduce_fill!(reduce, x, geoms, fill; kw...)
    else
        # Rasterize all features into a single bitarray.
        # The fill is the same so we can flatten the geometries together.
        bools = boolmask(dims(x, (XDim, YDim)), fc; kw...)
        # The fill `x` to match the masked values
        _fill!(x, bools, fill)
    end
    return x
end
function _rasterize!(x, ::GI.AbstractFeatureTrait, feature; fill, kw...)
    rasterize!(x, GI.geometry(feature); fill=_featurefillval(feature, fill), kw...)
end
function _rasterize!(x, ::GI.AbstractGeometryTrait, geom; fill, _buffer=nothing, kw...)
    ext = _extent(geom)
    x1 = view(x, Touches(ext))
    length(x1) > 0 || return x
    bools = if isnothing(_buffer)
        _init_bools(commondims(x1, (XDim, YDim)), Bool; missingval=false)
    else
        view(_buffer, Touches(ext))
    end
    boolmask!(bools, geom; kw...)
    _fill!(x1, bools, fill)
    return x
end
# Fill points
function _rasterize!(x, trait::GI.AbstractPointTrait, point; fill, kw...)
    _fill_point!(x, trait, point; fill, kw...)
    return x
end
# rasterize other iterables of features or gemoemtries
function _rasterize!(x, trait::Nothing, data; fill, reduce=last, kw...)
    # Check the itr length if we can, cycle if its 1
    fill_itr = if Base.IteratorSize(fill) isa Base.HasShape
        l = length(fill)
        if l == 1
            Iterators.cycle(fill)
        elseif l == n 
            fill
        else
            throw(ArgumentError("Length of fill $l does not match length of iterator $n"))
        end
    else
        fill
    end
    return _reduce_fill!(reduce, x, data, fill_itr; kw...)
end

# _reduce_fill!
#
# Mask `geoms` into each slice of a BitArray with the combined 
# dimensions of `x` and `geoms`, then apply the supplied reducing 
# function to each cell along the `:geoms` dimension. This uses a 
# conditional iterator over the fill values and the corresponding
# rasterized Bool - true included the fill value in the iterator, false excludes it.
#
# We get 64 Bool values to a regular `Int` meaning this doesn't scale too
# badly for large tables of geometries. 64k geometries and a 1000 * 1000
# raster needs 1GB of memory just for the `BitArray`.
function _reduce_fill!(f, x, geoms, fill_itr; kw...)
    # Define mask dimensions, the same size as the spatial dims of x
    spatialdims = commondims(x, (XDim, YDim))
    # Create a BitArray Raster with a dimension for the number of features to rasterize
    masks = boolmask(geoms; flat=false, kw...)
    spatialdims = commondims(x, (XDim, YDim))
    for ds in DimIndices(spatialdims)
        pixel_geom_list = view(masks, ds...)
        _apply_reduction!(f, x, ds, fill_itr, pixel_geom_list)
    end
    return x
end

# _apply_reduction!
#
# Apply a reducing functin over an iterable
# with performance optimisations where possible
function _apply_reduction!(f, x, ds, fill_itr, pixel_geom_list)
    any(pixel_geom_list) || return nothing
    # reduce is a custom reducing funcion
    iterator = (f for (f, b) in zip(fill_itr, pixel_geom_list) if b)
    pixel_value = f(iterator)
    x[ds...] = pixel_value
    return nothing
end
# last
function _apply_reduction!(f::typeof(last), x, ds, fill_itr, pixel_geom_list)
    local l = nothing
    for (v, b) in zip(fill_itr, pixel_geom_list)
        if b
            l = v
        end
    end
    if !isnothing(l)
        x[ds...] = l
    end
    return nothing
end
# first
function _apply_reduction!(f::typeof(first), x, ds, fill_itr, pixel_geom_list)
    for (v, b) in zip(fill_itr, pixel_geom_list)
        if b
            x[ds...] = v
            break
        end
    end
    return nothing
end
# count any all
for f in (:count, :any, :all)
    @eval function _apply_reduction!(f::typeof($f), x, ds, fill_itr, pixel_geom_list)
        x[ds...] = f(pixel_geom_list)
        return nothing
    end
end

_iscolumnfill(fill::Union{Symbol,Tuple{Symbol,Vararg}}) = true
_iscolumnfill(fill) = false

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
function _fill!(A::AbstractRaster{T}, B, fill, missingval=nothing) where T
    broadcast_dims!(A, A, B) do a, b
        val = b ? (fill isa Function ? fill(a) : fill) : a
        convert(T, val) # In case we are writing to disk
    end
end

function _at_or_contains(d, v, atol)
    selector = sampling(d) isa Intervals ? Contains(v) : At(v; atol=atol)
    DD.basetypeof(d)(selector)
end

_filter_name(name, fill::NamedTuple) = keys(fill)
_filter_name(name::NamedTuple, fill::NamedTuple) = keys(fill)
_filter_name(name::Nothing, fill::Nothing) = nothing
_filter_name(name::DimensionalData.NoName, fill::Union{Symbol,NTuple{<:Any,Symbol}}) = fill
_filter_name(name::Union{NamedTuple,Tuple,Array}, fill::NTuple{<:Any,Symbol}) = fill
function _filter_name(name::Union{NamedTuple,Tuple,Array}, fill::Union{Tuple,Array})
    length(name) == length(fill) || throw(ArgumentError("`name` keyword (possibly from `to` object) does not match length of fill. A fix is to use a `NamedTuple` for `fill`."))
    return name isa NamedTuple ? keys(name) : name
end
function _filter_name(name, fill)
    fill isa Union{Symbol,NTuple{<:Any,Symbol}} ? fill : name
end

function _dimkeys(table)
    Base.reduce(DEFAULT_TABLE_DIM_KEYS; init=()) do acc, key
        key in schema.names ? (acc..., key) : acc
    end
end
