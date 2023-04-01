function _auto_dim_columns(table, dims::Tuple)
    names = Tables.columnnames(table)
    if isempty(names)
        names = keys(first(Tables.rows(table)))
    end
    ds = DD.commondims(dims, (XDim, YDim, ZDim))
    return Tuple(DD.basetypeof(d)() for d in ds if DD.dim2key(d) in names)
end

_not_a_dimcol(table, dimcols::DimTuple) = _not_a_dimcol(table, map(DD.dim2key, dimcols))
_not_a_dimcol(table, dimcols::Tuple{Vararg{<:Pair}}) = _not_a_dimcol(table, map(last, dimcols))
_not_a_dimcol(table, dimcols::Tuple{}) = ()
function _not_a_dimcol(table, dimcols::Tuple{Vararg{Symbol}})
    dimcols = (dimcols..., :Band) # TODO deal with how annoying `Band` is 
    names = Tables.columnnames(table)
    if length(names) == 0
        names = keys(first(Tables.rows(table)))
    end
    not_dim_keys = Tuple(k for k in names if !(k in dimcols))
    return not_dim_keys
end
