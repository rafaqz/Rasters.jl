filter_ext(path, ext::AbstractString) = filter(fn -> splitext(fn)[2] == ext, readdir(path))
filter_ext(path, exts::Union{Tuple,AbstractArray}) = 
    filter(fn -> splitext(fn)[2] in exts, readdir(path))
filter_ext(path, ext::Nothing) = readdir(path)

# Check that array order matches expectation
checkarrayorder(A, order::Tuple) = map(checkarrayorder, dims(A), order)
checkarrayorder(dim::Dimension, order::Order) =
    arrayorder(dim) == order || @warn "Array order for `$(DD.basetypeof(order))` is `$(arrayorder(dim))`, usually `$order`"

checkorder(A, order::Tuple) = map(checkorder, dims(A), order)
checkorder(dim::Dimension, order::Order) =
    order(dim) == order || @warn "Array order for `$(DD.basetypeof(order))` is `$(order(dim))`, usually `$order`"

cleankeys(keys) = Tuple(map(Symbol, keys))

noindex_to_sampled(A) = rebuild(A; dims=noindex_to_sampled(dims(A)))
function noindex_to_sampled(dims::DimTuple)
    map(dims) do d
        lookup(d) isa NoLookup ? set(d, Sampled) : d
    end
end

function maybe_typemin_as_missingval(filename::String, A::AbstractRaster{T}) where T
    if ismissing(missingval(A))
        newmissingval = typemin(Missings.nonmissingtype(T)) 
        ext = splitext(filename)
        @warn "`missing` cant be written to $ext, typemin of $newmissingval used instead" 
        replace_missing(A, newmissingval)
    else
        A
    end
end

# We often need to convert the locus and the lookup in the same step,
# as doing it in the wrong order can give errors.
# function convert_locus_lookup(M1::Type{<:LookupArray}, L1::Type{<:Locus}, dim::Dimension)
#     _convert(S1, L1, sampling(dim), locus(dim), span(dim), dim)
# end

# _convert(::Type{M1}, ::Type{L1}, lookup::M2, l2::L2, span, dim) where {M1,M2<:M1,L1,L2<:L1} = dim
# _convert(::Type{M1}, ::Type{L1}, lookup::M2, l2::L2, span, dim) where {M1,M2<:M1,L1,L2} = shiftlocus(L1(), dim)
# _convert(::Type{M1}, ::Type{L1}, lookup::M2, l2::L2, span, dim) where {M1,M2<:M1,L1,L2} = 
#     _convert_by_locus(M1, L1, lookup, l2, span, dim)

# _convert_by_locus(M1, ::Type{Center}, lookup, l2::Union{Start,End}, span, dim) = 
#     _convert_by_lookup(M1, dim)
# _convert_by_locus(M1, L1::Type{Union{Start,End}}, lookup, l2::Center, span, dim) = 
#     _convert_by_lookup(M1, L1, dim)
# _convert_by_locus(M1, L1::Type{Start}, lookup, l2::End, span, dim) = 
#     convertlookup(M1, shiftlocus(L1, dim))
# _convert_by_locus(M1, L1::Type{End}, lookup, l2::Start, span, dim) = 
#     convertlookup(M1, shiftlocus(L1, dim))
# _convert_by_locus(M1, ::Type{L1}, lookup, l2::L2, span, dim) where {L1,L2<:L1} = 
#     convertlookup(M1, dim)

# # Projected will have an accurate center point on an equal scale, but Mapped may not. 
# # So we always shift the locus while in Projected lookup to avoid errors.
# _convert_by_lookup(::Type{Mapped}, dim) = convertlookup(Mapped, shiftlocus(Center(), dim))
# _convert_by_lookup(::Type{Projected}, dim) = shiftlocus(Center(), convertlookup(Projected, dim))
#

_not_a_dimcol(data, dimcols::DimTuple) = _not_a_dimcol(data, map(DD.dim2key, dimcols))
_not_a_dimcol(data, dimcols::Tuple{Vararg{<:Pair}}) = _not_a_dimcol(data, map(last, dimcols))
function _not_a_dimcol(data, dimcols::Tuple{Vararg{Symbol}})
    names = Tables.columnnames(data)
    if length(names) == 0
        names = keys(first(Tables.rows(data)))
    end
    not_dim_keys = Tuple(k for k in names if !(k in dimcols))
    return  not_dim_keys
end
