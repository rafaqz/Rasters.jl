
const BaseReduceFunc = Union{typeof(sum),typeof(prod),typeof(maximum),typeof(minimum),typeof(extrema)}

_reduce_op(::typeof(sum)) = Base.add_sum
_reduce_op(::typeof(prod)) = Base.mul_prod
_reduce_op(::typeof(minimum)) = min
_reduce_op(::typeof(maximum)) = max
_reduce_op(::typeof(extrema)) = Base._extrema_rf
_reduce_op(x) = nothing

_reduce_init(reduce, ::AbstractRaster{T}) where T = _reduce_init(reduce, T)
_reduce_init(reduce, st::AbstractRasterStack) = map(A -> _reduce_init(reduce, A), st)
_reduce_init(::Function, ::Type{T}) where T = zero(T)
_reduce_init(::typeof(sum), ::Type{T}) where T = zero(nonmissingtype(T))
_reduce_init(::typeof(prod), ::Type{T}) where T = oneunit(nonmissingtype(T))
_reduce_init(::typeof(minimum), ::Type{T}) where T = typemax(nonmissingtype(T))
_reduce_init(::typeof(maximum), ::Type{T}) where T = typemin(nonmissingtype(T))
_reduce_init(::typeof(extrema), ::Type{Union{Missing,T}}) where T<:Tuple{T1,T2} where {T1,T2} =
    (typemax(nonmissingtype(T1)), typemin(nonmissingtype(T2))) 
_reduce_init(::typeof(extrema), ::Type{T}) where T<:Tuple{T1,T2} where {T1,T2} =
    (typemax(nonmissingtype(T1)), typemin(nonmissingtype(T2))) 
_reduce_init(::typeof(extrema), ::Type{T}) where T<:Number =
    (typemax(nonmissingtype(T)), typemin(nonmissingtype(T)))

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
function rasterize(reduce::Function, data; kw...)
    rasterize(data; reduce, name=Symbol(string(reduce)), kw...)
end

const COUNT_NO_FILL = "`rasterize` with `count` does not use the `fill` keyword"
_count_fill(x) = x + 1
# Catch some functione early
#
# count is faster with an incrementing function as `fill`
function rasterize(reduce::typeof(count), data; fill=nothing, kw...)
    isnothing(fill) || @info COUNT_NO_FILL
    rasterize(data; name=:count, init=0, reduce=nothing, fill=_count_fill, kw...)
end
# `mean` is sum / count
# We can do better than this, but its easy for now
function rasterize(reduce::typeof(DD.Statistics.mean), data; fill, kw...)
    sums = rasterize(sum, data; fill, kw...)
    counts = rasterize(count, data; fill=nothing, kw...)
    rebuild(sums ./ counts; name=:mean)
end
# `any` is just boolmask
function rasterize(reduce::typeof(any), data; fill=nothing, kw...)
    isnothing(fill) || @info "`rasterize` with `any` does use the `fill` keyword"
    boolmask(data; kw..., name=:any, combine=true)
end
# `last` is the same as not reducing at all
function rasterize(reduce::typeof(last), data; fill, kw...)
    rasterize(data; name=:last, reduce=nothing, fill, kw...)
end

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
function _rasterize(to::DimTuple, data; fill, name=nothing, reduce=nothing, init=nothing, kw...)
    if Tables.istable(typeof(data)) # typeof so we dont check iterable table fallback in Tables.jl
        schema = Tables.schema(data)
        colnames = Tables.columnnames(Tables.columns(data))
        name = _filter_name(name, fill)
        # If fill is a symbol or tuple of Symbol we need to allocate based on the column type
        fillval = if fill isa Symbol
            fill in colnames || _fill_key_error(colnames, fill)
            zero(Tables.columntype(schema, fill))
        elseif fill isa Tuple{Symbol,Vararg}
            map(fill) do f
                f in colnames || _fill_key_error(colnames, fill)
                zero(Tables.columntype(schema, f))
            end
        elseif fill isa AbstractArray
            zero(eltype(fill))
        else
            # Otherwise allocate based on the fill type
            fill
        end
        init = if isnothing(init)
            isnothing(reduce) ? nothing : _reduce_init(reduce, typeof(fillval))
        else
            init
        end
        return _create_rasterize_dest(to; fill=fillval, name, init, kw...) do dest
            rasterize!(dest, data; reduce, fill, init, kw...)
        end
    else
        _rasterize(to, GeoInterface.trait(data), data; reduce, fill, init, name, kw...)
    end
end
function _rasterize(to::DimTuple, ::GI.AbstractFeatureCollectionTrait, fc; name, fill, kw...)
    # TODO: how to handle when there are fillvals with different types
    fillval = _featurefillval(GI.getfeature(fc, 1), fill)
    name = _filter_name(name, fill)
    return _create_rasterize_dest(to1; fill=fillval, name, kw...) do dest
        rasterize!(dest, fc; fill, kw...)
    end
end
function _rasterize(to::DimTuple, ::GI.AbstractFeatureTrait, feature; fill, name, kw...)
    fillval = _featurefillval(feature, fill)
    name = _filter_name(name, fill)

    return _create_rasterize_dest(to; fill=fillval, name, kw...) do dest
        rasterize!(dest, feature; fill, kw...)
    end
end
function _rasterize(to::DimTuple, ::Union{Nothing,GI.AbstractGeometryTrait}, data; fill, kw...)
    return _create_rasterize_dest(to; fill, kw...) do dest
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
function rasterize!(reduce::typeof(count), x::RasterStackOrArray, data; fill=nothing, kw...)
    isnothing(fill) || @info COUNT_NO_FILL
    rasterize!(x::RasterStackOrArray, data; kw..., reduce=nothing, op=nothing, fill=_count_fill)
end
function rasterize!(x::RasterStackOrArray, data; fill, reduce=nothing, kw...)
    function _rasterize_point_table_inner!(x, pointcols, fill_itr; kw...)
        for (point, fill) in  zip(zip(pointcols...), fill_itr)
            _fill_point!(x, point; fill, kw...)
        end
    end
    # Check if it's a Tables.jl compatible object
    if Tables.istable(typeof(data))
        schema = Tables.schema(data)
        geomcolname = first(GI.geometrycolumns(data))
        cols = Tables.columns(data)

        if geomcolname in Tables.columnnames(cols)
            # Its a geometry table
            geomcol = Tables.getcolumn(cols, geomcolname)
            fill_itr = _iterable_fill(cols, fill)
            _rasterize!(x, nothing, geomcol; reduce, fill=fill_itr, kw...)
        else
            # Its a point table
            dimscols = _auto_dim_columns(data, dims(x))
            pointcols = map(k -> Tables.getcolumn(data, k), map(DD.dim2key, dimscols))
            reduce == last || throw(ArgumentError("Can only reduce with `last` on point tables. Make a github issue at the Rasters.jl repository if you need something else."))
            _buffer = _init_bools(x, Bool, geom; kw..., missingval=false)
            _rasterize_point_table_inner!(x, pointcols, fill_itr; _buffer, reduce, kw...)
        end
    else
        # Otherwise maybe a Geometry, FeatureCollection or some iterator
        _rasterize!(x, GI.trait(data), data; fill, reduce, kw...)
    end
end

function _rasterize!(x, ::GI.AbstractFeatureCollectionTrait, fc; reduce=nothing, fill, kw...)
    fill_itr = _iterable_fill(fc, keyorfill)
    features = GI.feaatures(fc)
    _rasterize!(x, nothing, features; fill=fill_itr, kw...)
    return x
end

# We rasterize all iterables here
function _rasterize!(x, trait::Nothing, geoms::AbstractVector;
    fill, reduce=nothing, op=nothing, _buffer=nothing, kw...
)
    # Check the itr length if we can, cycle if its 1
    fill_itr = _iterable_fill(geoms, fill)

    # Special-casing for performance
    if isnothing(reduce) && isnothing(op)
        if fill_itr isa Iterators.Cycle
            # We dont need to iterate the fill, so just mask
            mask = boolmask(geoms; to=x, combine=true, kw...)
            # And broadcast the fill
            broadcast_dims!(x, x, mask) do v, m
                m ? fill_itr.xs : v
            end
        elseif fill_itr isa NamedTuple
            if all(f -> f isa Iterators.Cycle, fill_itr)
                # We dont need to iterate the fill, so just mask
                mask = boolmask(geoms; to=x, combine=true, kw...)
                foreach(x, fill_itr) do A, f 
                    # And broadcast the fill
                    broadcast_dims!(A, A, mask) do v, m
                        m ? f.xs : v
                    end
                end
            else
                Threads.@threads for i in eachindex(geoms)
                    geom = geoms[i]
                    fill = NamedTuple{keys(fill_itr)}(map(f -> _getfill(f, i), fill_itr))
                    _rasterize!(x, GI.trait(geom), geom; fill, kw...)
                end
            end
        else
            Threads.@threads for i in eachindex(geoms)
                geom = geoms[i]
                fill = _getfill(fill_itr, i)
                _rasterize!(x, GI.trait(geom), geom; kw..., fill)
            end
        end
    else 
        # See if there is a reducing operator passed in, or matching `reduce`
        op = isnothing(op) ? _reduce_op(reduce) : op
        init = _reduce_init(reduce, x)
        if isnothing(op)
            # If there still isn't any `op`, reduce over all the values later rather than iteratively
            _reduce_fill!(reduce, x, geoms, fill_itr; kw..., init)
        else # But if we can, use op in a fast iterative reduction
            Threads.@threads for i in eachindex(geoms)
                geom = geoms[i]
                fill = _getfill(fill_itr, i)
                _rasterize!(x, GI.trait(geom), geom; kw..., fill, op, init)
            end
        end
    end
    return x
end

function _rasterize!(x, ::GI.AbstractFeatureTrait, feature; fill, kw...)
    rasterize!(x, GI.geometry(feature); fill=_featurefillval(feature, fill), kw...)
end
function _rasterize!(x, ::GI.AbstractGeometryTrait, geom; fill, _buffer=nothing, op=nothing, init=nothing, missingval=missing, kw...)
    ext = _extent(geom)
    x1 = view(x, Touches(ext))
    length(x1) > 0 || return x

    bools = if isnothing(_buffer)
        _init_bools(x1, Bool, geom; kw..., missingval=false)
    else
        b = view(_buffer, Touches(ext))
    end
    boolmask!(bools, geom; kw...)
    _fill!(x1, bools, fill, op, init, missingval)
    return x
end
# Fill points
function _rasterize!(x, trait::GI.AbstractPointTrait, point; fill, kw...)
    _fill_point!(x, trait, point; fill, kw...)
    return x
end

# _reduce_fill!
#
# Mask `geoms` into each slice of a BitArray with the combined
# dimensions of `x` and `geoms`, then apply the supplied reducing
# function to each cell along the `:geometry` dimension. This uses a
# conditional iterator over the fill values and the corresponding
# rasterized Bool - true included the fill value in the iterator, false excludes it.
#
# We get 64 Bool values to a regular `Int` meaning this doesn't scale too
# badly for large tables of geometries. 64k geometries and a 1000 * 1000
# raster needs 1GB of memory just for the `BitArray`.
function _reduce_fill!(f, st::AbstractRasterStack, geoms, fill_itr::NamedTuple{K}; kw...) where K
    # Define mask dimensions, the same size as the spatial dims of x
    spatialdims = commondims(st, DEFAULT_POINT_ORDER)
    # Mask geoms as separate bool layers
    masks = boolmask(geoms; to=st, combine=false, kw...)
    # Use a generator over the array axis in case the iterator has no length
    geom_axis = axes(masks, Dim{:geometry}())
    fill = map(itr -> [v for (_, v) in zip(geom_axis, itr)], fill_itr)
    _reduce_fill_inner!(f, st, geoms, fill, masks, parent(parent(masks)))
end
function _reduce_fill!(f, A::AbstractRaster, geoms, fill_itr; kw...)
    # Define mask dimensions, the same size as the spatial dims of x
    spatialdims = commondims(A, DEFAULT_POINT_ORDER)
    # Mask geoms as separate bool layers
    masks = boolmask(geoms; to=A, combine=false, kw...)
    # Use a generator over the array axis in case the iterator has no length
    geom_axis = parent(axes(masks, Dim{:geometry}()))
    fill = [val for (i, val) in zip(geom_axis, fill_itr)]
    _reduce_fill_inner!(f, A, geoms, fill, masks)
end

# Separated as a function barrier for type stability
function _reduce_fill_inner!(f, x::AbstractRaster{T}, geoms, fill::Union{AbstractVector,NamedTuple}, masks) where T
    Threads.@threads for D in DimIndices(dims(x, DEFAULT_POINT_ORDER))
        # Do DimensionalData.jl indexing manually to avoid taking a view of the index
        I = dims2indices(masks, D)
        newval = _apply_reduction!(f, fill, view(parent(parent(masks)), I...))::Union{T,Nothing}
        if !isnothing(newval)
            @inbounds x[D...] = newval::T
        end
    end
    return x
end

# _apply_reduction!
#
# Apply a reducing functin over an iterable
# with performance optimisations where possible
@inline function _apply_reduction!(f, fill_itr, pixel_geom_list)
    any(pixel_geom_list) || return nothing
    iterator = (f for (f, b) in zip(fill_itr, pixel_geom_list) if b)
    return f(iterator)
end

# _featurefillval
# Get fill value from a feature, or use fill itself
_featurefillval(feature, fill::Nothing) = first(GI.properties(feature))
_featurefillval(feature, fill::Symbol) = GI.properties(feature)[fill]
_featurefillval(feature, fill::Val) = _featurefillval(feature, _unwrap(fill))
_featurefillval(feature, fill::NamedTuple) = map(f -> _featurefillval(feature, f), _unwrap(fill))
function _featurefillval(feature, fill::NTuple{<:Any,Symbol})
    map(fill) do key
        getproperty(GI.properties(feature), key)
    end
end
_featurefillval(feature, fill) = fill

function _create_rasterize_dest(f, dims; fill, name=nothing, init=nothing, kw...)
    _create_rasterize_dest(f, fill, init, name, dims; fill, kw...)
end
function _create_rasterize_dest(f, fill::Union{Tuple,NamedTuple}, init, keys, dims; kw...)
    _create_rasterize_dest(f, fill, init, DD.uniquekeys(fill), dims; kw...)
end
function _create_rasterize_dest(f, fill::Union{Tuple,NamedTuple}, init, keys::Union{Tuple,NamedTuple}, dims;
    filename=nothing, missingval=nothing, metadata=NoMetadata(), suffix=nothing, kw...
)
    dims = _as_intervals(dims) # Only makes sense to rasterize to intervals
    init1 = isnothing(init) ? map(_ -> nothing, fill) : init
    layers = map(keys, values(fill), values(init1)) do key, fillval, initval
        missingval = isnothing(missingval) ? _writeable_missing(filename, typeof(fillval)) : missingval
        _alloc_rasterize(filename, fillval, initval, dims; name, metadata, missingval, suffix=key) do a
            a
        end
    end
    st = RasterStack(layers, dims; keys, metadata)
    open(f, RasterStack(layers, dims; keys, metadata))
    return st
end
function _create_rasterize_dest(f, fill, init, name, dims;
    filename=nothing, missingval=nothing, metadata=NoMetadata(), suffix=nothing, kw...
)
    dims = _as_intervals(dims) # Only makes sense to rasterize to intervals
    missingval = isnothing(missingval) ? _writeable_missing(filename, typeof(val)) : missingval
    A = _alloc_rasterize(filename, fill, init, dims; name, metadata, missingval, suffix) do a
        f(a)
    end
    return Raster(A, dims; name, missingval, metadata)
end

function _alloc_rasterize(f, filename, fill, init, to; missingval, kw...)
    T = if isnothing(init)
        if fill isa AbstractArray
            eltype(fill)
        elseif fill isa Function
            typeof(fill(missingval))
        else
            typeof(fill)
        end
    else
        typeof(init)
    end
    T1 = promote_type(typeof(missingval), T)
    A = create(filename, T1, to; missingval, kw...)
    # TODO f should apply to the file when it is initially created
    # instead of reopening but we need a `create(f, filename, ...)` method
    open(A; write=true) do A
        A .= Ref(missingval)
        f(A)
    end
    return A
end


function _fill!(st::AbstractRasterStack, B, fill, op, init, missingval)
    if fill isa NamedTuple
        init = isnothing(init) ? map(_ -> nothing, fill) : init
        foreach(values(st), fill, init) do a, f, i
            _fill!(a, B, f, op, i, missingval)
        end
    end
    return st
end
# If the array is initialised, we can use the existing values
function _fill!(A::AbstractRaster{T}, B, fill, op, init, missingval) where T
    broadcast_dims!(A, A, B) do a, b
        convert(T, _choose_fill(a, b, fill, op, init, missingval))
    end
end

function _choose_fill(a, b, fill, op, init, missingval)

    if b
        newval = if isnothing(op) && !(fill isa Function)
            fill
        else
            # Either a reducing `op` or a function of the current value
            # Get the value: a, init, fill or fill(a)
            a1 = if ismissing(a)
                if ismissing(missingval) & !isnothing(init)
                    init
                else
                    a
                end
            elseif ismissing(missingval)
                fill isa Function ? fill(a) : a
            elseif a == missingval
                fill isa Function ? fill(init) : init
            else
                fill isa Function ? fill(a) : a
            end

            if isnothing(op)
                a1
            else
                if op === Base._extrema_rf
                    op(a1, (fill, fill))
                else
                    op(a1, fill)
                end
            end
        end
    else
        a
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

# A Tuple of `Symbol` is multiple keys to make a RasterStack
_iterable_fill(data, keys::Tuple{Symbol,Vararg}) =
    NamedTuple{keys}(map(k -> _iterable_fill(data, k), keys))
# A Symbol is a Table or FeatureCollection key, it cant be used as fill itself
function _iterable_fill(data, key::Symbol)
    # For column tables, get the column now
    if !Tables.isrowtable(typeof(data))
        names = Tables.columnnames(Tables.columns(data))
        key in names || _fill_key_error(names, key)
        Tables.getcolumn(data, key)
    else
        # For row tables and FeatureCollection, get the values row by row
        # We lift the Symbols to the type domain here so the generator is type stable later
        Val{key}()
    end
end
_iterable_fill(data, fill::NamedTuple) = map(f -> _iterable_fill(data, f), fill)
# Inspect our data and fill as much as possible to check they match
# and cycle any fill of known length one
function _iterable_fill(data, fill)
    if Base.IteratorSize(fill) isa Base.HasShape
        l = length(fill)
        if l == 1
            # Cycle all length one iterables to fill every row
            Iterators.cycle(fill)
        elseif Base.IteratorSize(data) isa Union{Base.HasShape,Base.HasLength}
            # We know the data and iterator length, so check that they match to catch errors early with a clean message
            n = length(data)
            if l == n
                fill
            else
                throw(ArgumentError("Length of fill $l does not match length of iterator $n"))
            end
        else
            # We don't know the length of the data so let it error later if it has to
            fill
        end
    else
        # We don't knwo the length of the data or the fill, so whatever happens, happens
        fill
    end
end

_getfill(itrs::NamedTuple, i::Int) = map(itr -> _getfill(itr, i), itrs)
_getfill(itr::AbstractArray, i::Int) = itr[i]
_getfill(itr::Iterators.Cycle, i::Int) = first(itr)
_getfill(itr, i) = itr
