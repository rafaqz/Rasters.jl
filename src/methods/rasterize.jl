struct _Undefined end
struct _Defined end

"""
    rasterize(x; kw...)
    rasterize(points, values; kw...)

Rasterize the points and values in a Tables.jl compatible object,
or separate `points` and `values` vetors/iterators, which must be
the same length.

If a GeoInterface `AbstractGeometry` or nested `Vector`s of `Tuple`/`Vector`
points is passed in, a `fill` keyword is also required to provide the value
that the array will be filled with. If `fill` is a function, it will be
applied to the existing value present in the array.

# Arguments

- `x`: a Tables.jl compatible object containing points and values columns,
    an GeoInterface.jl `AbstractGeometry`, or a nested `Vector` of `Vectors`.
- `points`: A `Vector` or nested `Vectors` holding `Vector` or `Tuple` of `Real`
- `values` A `Vector` of values to be written to a `Raster`, or a Vector of `NamedTuple`
    to write to a `RasterStack`.

# Keywords

These are detected automatically from `A` and `data` where possible.

- `to`: a `Raster`, `RasterStack` of `Tuple` of `Dimension` to use as a to.
- `order`: A `Tuple` of pairs `Dim => Symbol` for the keys in the data that match
    the dimension, or for the order of dimensions in ppoint values, like `(X, Y)`.
- `atol`: an absolute tolerance for rasterizing to dimensions with `Points` sampling.
- `filename`: a filename to write to directly, useful for large files.
- `suffix`: a string or value to append to the filename.
    A tuple of `suffix` will be applied to stack layers. `keys(st)` are the default.

## Geometry keywords

These can be used when a `GeoInterface.AbstractGeometry` is passed in.

- `fill`: the value to fill a polygon with, if `data` is a polygon. 
- `shape`: Force `data` to be treated as `:polygon`, `:line` or `:point`.

And specifically for `shape=:polygon`:

- `boundary`: include pixels where the `:center` is inside the polygon, where 
    the line `:touches` the pixel, or that are completely `:inside` inside the polygon.

## Table keywords

- `name`: A `Symbol` to return a `Raster` from a single column,
    or `Tuple` of `Symbol` to return a `RasterStack` from multiple columns.

# Example

Rasterize a shapefile for China and plot, with a border.

```jldoctest
using Rasters, Plots, Dates, Shapefile, Downloads
using Rasters.LookupArrays

# Download a borders shapefile
shapefile_url = "https://github.com/nvkelso/natural-earth-vector/raw/master/10m_cultural/ne_10m_admin_0_countries.shp"
shapefile_name = "boundary_lines.shp"
isfile(shapefile_name) || Downloads.download(shapefile_url, shapefile_name)

# Loade the shapes for china
china_border = Shapefile.Handle(shapefile_name).shapes[10]

# Define dims for the china area
dms = Y(Projected(15.0:0.1:55.0; order=ForwardOrdered(), span=Regular(0.1), sampling=Intervals(Start()), crs=EPSG(4326))), 
      X(Projected(70.0:0.1:140; order=ForwardOrdered(), span=Regular(0.1), sampling=Intervals(Start()), crs=EPSG(4326)))

# Rasterize the border polygon 
using BenchmarkTools
using ProfileView
@time china = rasterize(china_border; 
    to=dms, missingval=-9999, fill=1, 
    order=(X, Y),# shape=:line, 
    boundary=:touches,
    # filename="rasterize.tif"
)

# And plot
p = plot(china; color=:spring)
plot!(p, china_border; fillalpha=0, linewidth=0.6)
savefig("build/china_rasterized.png")

# output
```

![rasterize](china_rasterized.pngfill)

$EXPERIMENTAL
"""
function rasterize(data; to, order=nothing, name=nothing, kw...) 
    if Tables.istable(data)
        points, values, order, name = _table_to_points_values(to, data; order, name)
        return _rasterize(to, points, values; order, name, kw...)
    else
        order = isnothing(order) ? DEFAULT_POINT_ORDER : order
        order = dims(to, order)
        return _rasterize(to, data; order, init=_Undefined(), kw...)
    end
end
function rasterize(points, values; to, kw...) 
    return _rasterize(to, points, values; kw...)
end

function _rasterize(to::DimTuple, points;
    filename=nothing, suffix=nothing, metadata=NoMetadata(), name=nothing, parent=nothing,
    fill, eltype=typeof(fill), missingval=_writeable_missing(filename, eltype), kw...
)
    A = _alloc_rasterize(filename, eltype, to; missingval, suffix, parent) do a
        rasterize!(a, points; fill, missingval, kw...)
    end
    return A
end
function _rasterize(to::DimTuple, points, vals;
    missingval=nothing, filename=nothing, suffix=nothing,
    metadata=NoMetadata(), name=nothing, parent=nothing, kw...
)
    firstval = first(vals)
    if firstval isa NamedTuple
        # Rasterize mutiple values to a stack
        layers = map(keys(firstval), values(firstval)) do key, val
            missingval = isnothing(missingval) ? _writeable_missing(filename, typeof(val)) : missingval
            _alloc_rasterize(filename, typeof(val), to; name, metadata, missingval, suffix=key, parent) do a
                a .= missingval
            end
        end |> NamedTuple{keys(firstval)}
        st = RasterStack(layers, to)
        return rasterize!(st, points, vals; kw...)
    end
    # Rasterize to an array
    missingval = isnothing(missingval) ? _writeable_missing(filename, typeof(firstval)) : missingval
    A = _alloc_rasterize(filename, typeof(firstval), to; name, metadata, missingval, suffix, parent) do a
        a .= missingval
    end
    return rasterize!(A, points, vals; kw...)
end
function _rasterize(to::AbstractRaster, args...; missingval=missingval(to), kw...)
    _rasterize(dims(to), args...; missingval, kw...)
end
function _rasterize(to::AbstractRasterStack, args...; kw...)
    return _rasterize(dims(to), args...; kw...)
end

function _alloc_rasterize(f, filename, T, to; missingval, suffix=nothing, kw...)
    T = promote_type(typeof(missingval), T)
    A = create(filename, T, to; suffix, missingval, kw...)
    return A
end

"""
    rasterize!(x, data; order, name, atol)
    rasterize!(x, points, vals; order, atol)

Rasterize the points and vals in `data`, or the `points` and `vals` objects,
into the [`Raster`](@ref) or [`RasterStack`](@ref) `x`.

# Arguments

- `x`: a `Raster` or `RasterStack` to rasterize to.
- `data`: a Tables.jl compatible object containing points and values or a
    polygon - an GeoInterface.jl `AbstractGeometry`, or a nested `Vector` of `Vectors`.
- `points`: A `Vector` or nested `Vector` holding `Vector` or `Tuple` of `Real`
- `vals` A `Vector` of values to be written when `x` is a `Raster`, or a Vector of
    `NamedTuple` to write when `x` is a `RasterStack`.

# Keywords

These are detected automatically from `A` and `data` where possible.

- `order`: A `Tuple` of pairs `Dim => Symbol` for the keys in the data that match
    the dimension, or for the order of dimensions in points, like `(X, Y)`.
- `atol`: an absolute tolerance for rasterizing to dimensions with `Points` sampling.
- `filename`: a filename to write to directly, useful for large files.
- `suffix`: a string or value to append to the filename.
    A tuple of `suffix` will be applied to stack layers. `keys(st)` are the default.

## Geometry keywords

These can be used when a `GeoInterface.AbstractGeometry` is passed in.

- `fill`: the value to fill a polygon with, if `data` is a polygon.
    `fill can also be a `Function` of the existing value.
- `shape`: Force `data` to be treated as `:polygon`, `:line` or `:point`.

And specifically for `shape=:polygon`:

- `boundary`: include pixels where the `:center` is inside the polygon, where 
    the line `:touches` the pixel, or that are completely `:inside` inside the polygon.

## Table keywords

- `name`: A `Symbol` to return a `Raster` from a single column,
    or `Tuple` of `Symbol` to return a `RasterStack` from multiple columns.

# Example

```jldoctest
using Rasters, Plots, Dates, Shapefile, GeoInterface, Downloads
using Rasters.LookupArrays

# Download a borders shapefile
shapefile_url = "https://github.com/nvkelso/natural-earth-vector/raw/master/10m_cultural/ne_10m_admin_0_countries.shp"
shapefile_name = "boundary_lines.shp"
isfile(shapefile_name) || Downloads.download(shapefile_url, shapefile_name)

# Load the shapes for denmark
indonesia_border = Shapefile.Handle(shapefile_name).shapes[1]

# Make an empty EPSG 4326 projected Raster of the area of Indonesia
dimz = Y(-15.0:0.1:10.9; mode=Projected(; sampling=Intervals(Start()), crs=EPSG(4326))), 
       X(90.0:0.1:145; mode=Projected(; sampling=Intervals(Start()), crs=EPSG(4326)))

A = Raster(zeros(UInt16, dimz); missingval=0)

# Rasterize each island with a different number
using ProfileView
@profview for (i, shp) in enumerate(coordinates(indonesia_border))
    rasterize!(A, shp; fill=i, order=(X, Y), shape=:polygon, boundary=:center)
end
@time for (i, shp) in enumerate(coordinates(indonesia_border))
    rasterize!(A, shp; fill=i, order=(X, Y), shape=:polygon, boundary=:all)
end

# And plot
p = plot(A; color=:spring)
plot!(p, indonesia_border; fillalpha=0, linewidth=0.7)
savefig("build/indonesia_rasterized.png")

# output

```

![rasterize](indonesia_rasterized.png)

$EXPERIMENTAL
"""
function rasterize!(x::RasterStackOrArray, data;
    order=nothing, name=nothing, missingval=missingval(x), kw...
)
    if Tables.istable(data)
        _rasterize_table!(x, data; name, kw...)
    else
        order = isnothing(order) ? DEFAULT_POINT_ORDER : order
        buffer = Raster(falses(commondims(x, order)))
        _rasterize_geometry!(x, data; order, buffer, init=_Defined(), missingval, kw...)
    end
end
function rasterize!(A::AbstractRaster, geom::GI.AbstractGeometry, vals; shape=:point, kw...)
    rasterize!(A, _flat_nodes(geom), vals; kw...)
end
function rasterize!(st::AbstractRasterStack, geom::GI.AbstractGeometry, vals; shape=:point, kw...)
    rasterize!(st, _flat_nodes(geom), vals; kw...)
end
function rasterize!(x::AbstractRasterStack, points, vals; name=nothing, kw...)
    firstval = first(vals)
    name = if isnothing(name) 
        if firstval isa NamedTuple
            keys(first(vals)) 
        else
            length(name) == keys(x) || throw(ArgumentError("`name` keyword does not match number of layers in stack"))
            keys(x)
        end
    else
        if firstval isa NamedTuple
            name == keys(first(vals)) || throw(ArgumentError("`name` keyword does not match point names"))
        elseif firstval isa Union{Tuple,Array}
            length(name) == length(firstval) || throw(ArgumentError("`name` keyword does not match length of point"))
        end
        name
    end
    return _rasterize!(x[name], points, vals; name, kw...)
end
function rasterize!(A::AbstractRaster, points, vals; name=nothing, kw...)
    _rasterize!(A, points, vals; name, kw...)
end

function _rasterize!(x::RasterStackOrArray, points, vals; 
    order=DEFAULT_POINT_ORDER, atol=nothing, name=nothing
)
    check_points(points, order)
    isdisk(first(x)) && _warn_disk()
    ordered_dims = dims(x, order)
    _without_mapped_crs(x) do x
        map(points, vals) do p, v
            any(map(ismissing, p)) && return nothing
            selectors = map((d, x) -> _at_or_contains(d, x, atol), ordered_dims, p)
            all(map(DD.hasselection, ordered_dims, selectors)) || return nothing
            _set_at_selection!(x, selectors, v)
            return nothing
        end
    end
    return x
end

function _set_at_selection!(st::AbstractRasterStack, selectors, vals)
    map(values(st), values(vals)) do A, val
        _set_at_selection!(A, selectors, val)
    end
end
function _set_at_selection!(A::AbstractRaster, selectors, val::Union{NamedTuple,Tuple})
    _set_at_selection!(A, selectors, first(val))
end
function _set_at_selection!(A::AbstractRaster, selectors, val)
    if length(selectors) == length(dims(A))
        A[selectors...] = val
    else
        A[selectors...] .= val
    end
end

function check_points(points, order)
    lp = length(first(points))
    lo = length(order) 
    if lp !== lo 
        throw(ArgumentError("length of points in `points` ($lp) is different to the length of `order` ($lo)")) 
    end
end

# _rasterize_table!
function _rasterize_table!(x::RasterStackOrArray, data; order=nothing, name=nothing, kw...)
    points, vals, order, name = _table_to_points_values(x, data; order, name)
    return rasterize!(x, points, vals; order, name, kw...)
end

# _rasterize_geometry!
function _rasterize_geometry!(x::RasterStackOrArray, geom; 
    buffer, fill, order, init, missingval, kw...
)
    boolmask!(buffer, geom; order, kw...)
    return _fill!(x, buffer, fill, init, missingval)
end

_fill(B, fill, missingval) where T = broadcast_dims(x -> x ? fill : missingval, B)

function _fill!(A::AbstractRasterStack, B, fill, args...)
    map((a, f) -> _fill!(a, B, f, args...), st, fill)
    return A
end
# If the array is initialised, we can use the existing values
function _fill!(A::AbstractRaster{T}, B, fill, init::_Defined, missingval) where T
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
