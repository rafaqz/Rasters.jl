function _table_to_points_values(to, st::RasterStack; name=nothing, order=nothing, kw...)
    name = isnothing(name) ? keys(st) : name
    order = isnothing(order) ? map(d -> basetypeof(d) => dim2key(d), dims(st)) : order
    _table_to_points_values(to, DimTable(st); name, order, kw...)
end
function _table_to_points_values(to, A::Raster; name=nothing, order=nothing, kw...)
    name = isnothing(name) ? DD.name(A) : name
    order = isnothing(order) ? map(d -> basetypeof(d) => dim2key(d), dims(A)) : order
    _table_to_points_values(to, DimTable(A); name, order, kw...)
end
function _table_to_points_values(to, table; order, name)
    order = _table_point_order(dims(to), table, order)
    points = _table_points(table, order)
    dimcols = map(last, order)
    name = isnothing(name) ? _not_a_dimcol(table, order) : name
    vals = if name isa Symbol
        row_iter(Val{name}(), table)
    elseif name isa NTuple{<:Any,Symbol}
        row_iter(NamedTuple{name}(name), table)
    else
        throw(ArgumentError("`name` must be a `Symbol` or `Tuple` of `Symbol`. Got $name"))
    end
    dimorder = map(first, order)
    return points, vals, dimorder, name
end

@noinline row_iter(::Val{N}, table) where N = (r[N] for r in Tables.rows(table))
@noinline function row_iter(::NamedTuple{N}, table) where N
    (NamedTuple{N}(map(n -> r[n], N)) for r in Tables.rows(table))
end

function _table_point_order(dims::DimTuple, table, order::Tuple{<:Pair,Vararg{<:Pair}})
    all(hasdim(dims, map(first, order))) || throw(ArgumentError("Order dims $(otherdims(dims, order)) not found in dims"))
    names = Tables.columnnames(table)
    if names === ()
        names = keys(first(Tables.rows(table)))
    end
    foreach(order) do p
        last(p) in names || throw(ArgumentError("$(last(order)) not found in table columns $names"))
    end
    return order
end
function _table_point_order(dims::DimTuple, table, order)
    msg = "Specify `order` as `dim => column` pairs, e.g. `order=(X => :xcol, Y => :ycol)`. Got $order"
    throw(ArgumentError(msg))
end
function _table_point_order(dims::DimTuple, table, order::Nothing)
    order = _auto_dim_columns(dims, table)
    names = Tables.columnnames(table)
    if names === ()
        names = keys(first(Tables.rows(table)))
    end
    if order === ()
        msg = "Could not detect dimension columns in `$names`. Please specify `order` as `dim => :column` pairs, e.g. `order=(X => :xcol, Y => :ycol)`"
        throw(ArgumentError(msg))
    end
    return order
end

function _table_points(table, order)
    isdisk(table) && _warn_disk()
    dimnames = map(last, order)
    points = (map(n -> r[n], dimnames) for r in Tables.rows(table))
    return points
end

function _table_bounds(table, order)
    points = _table_points(table, order)
    bounds = map(ntuple(x -> x, length(order))) do i
        extrema(x -> x[i], points)
    end
    return bounds
end

function _wrapped_table_bounds(dims::DimTuple, table, order::Tuple)
    ds = DD.dims(dims,  map(first, order))
    map(rebuild, ds, _table_bounds(table, order))
end

function _auto_dim_columns(dims::DimTuple, table)
    names = Tables.columnnames(table)
    if isempty(names)
        names = keys(first(Tables.rows(table)))
    end
    ds = DD.commondims(dims, (XDim, YDim, ZDim))
    return Tuple(DD.basetypeof(d)() => DD.dim2key(d) for d in ds if DD.dim2key(d) in names)
end

_not_a_dimcol(data, dimcols::DimTuple) = _not_a_dimcol(data, map(DD.dim2key, dimcols))
_not_a_dimcol(data, dimcols::Tuple{Vararg{<:Pair}}) = _not_a_dimcol(data, map(last, dimcols))
_not_a_dimcol(data, dimcols::Tuple{}) = ()
function _not_a_dimcol(data, dimcols::Tuple{Vararg{Symbol}})
    dimcols = (dimcols..., :Band) 
    names = Tables.columnnames(data)
    if length(names) == 0
        names = keys(first(Tables.rows(data)))
    end
    not_dim_keys = Tuple(k for k in names if !(k in dimcols))
    return not_dim_keys
end

