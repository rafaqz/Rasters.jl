filter_ext(path, ext::AbstractString) = filter(fn -> splitext(fn)[2] == ext, readdir(path))
filter_ext(path, exts::Union{Tuple,AbstractArray}) = 
    filter(fn -> splitext(fn)[2] in exts, readdir(path))
filter_ext(path, ext::Nothing) = readdir(path)

# Check that array order matches expectation
checkarrayorder(A, order::Order) = map(d -> checkarrayorder(d, order), dims(A))
checkarrayorder(A, order::Tuple) = map(checkarrayorder, dims(A), order)
checkarrayorder(dim::Dimension, order::Order) =
    arrayorder(dim) == order || @warn "Array order for `$(DD.basetypeof(order))` is `$(arrayorder(dim))`, usually `$order`"

checkindexorder(A, order::Order) = map(d -> checkindexorder(d, order), dims(A))
checkindexorder(A, order::Tuple) = map(checkindexorder, dims(A), order)
checkindexorder(dim::Dimension, order::Order) =
    indexorder(dim) == order || @warn "Array order for `$(DD.basetypeof(order))` is `$(indexorder(dim))`, usually `$order`"

cleankeys(keys) = Tuple(map(Symbol, keys))

# We often need to convert the locus and the mode in the same step,
# as doing it in the wrong order can give errors.
# function convert_locus_mode(M1::Type{<:IndexMode}, L1::Type{<:Locus}, dim::Dimension)
#     _convert(S1, L1, sampling(dim), locus(dim), span(dim), dim)
# end

# _convert(::Type{M1}, ::Type{L1}, mode::M2, l2::L2, span, dim) where {M1,M2<:M1,L1,L2<:L1} = dim
# _convert(::Type{M1}, ::Type{L1}, mode::M2, l2::L2, span, dim) where {M1,M2<:M1,L1,L2} = shiftlocus(L1(), dim)
# _convert(::Type{M1}, ::Type{L1}, mode::M2, l2::L2, span, dim) where {M1,M2<:M1,L1,L2} = 
#     _convert_by_locus(M1, L1, mode, l2, span, dim)

# _convert_by_locus(M1, ::Type{Center}, mode, l2::Union{Start,End}, span, dim) = 
#     _convert_by_mode(M1, dim)
# _convert_by_locus(M1, L1::Type{Union{Start,End}}, mode, l2::Center, span, dim) = 
#     _convert_by_mode(M1, L1, dim)
# _convert_by_locus(M1, L1::Type{Start}, mode, l2::End, span, dim) = 
#     convertmode(M1, shiftlocus(L1, dim))
# _convert_by_locus(M1, L1::Type{End}, mode, l2::Start, span, dim) = 
#     convertmode(M1, shiftlocus(L1, dim))
# _convert_by_locus(M1, ::Type{L1}, mode, l2::L2, span, dim) where {L1,L2<:L1} = 
#     convertmode(M1, dim)

# # Projected will have an accurate center point on an equal scale, but Mapped may not. 
# # So we always shift the locus while in Projected mode to avoid errors.
# _convert_by_mode(::Type{Mapped}, dim) = convertmode(Mapped, shiftlocus(Center(), dim))
# _convert_by_mode(::Type{Projected}, dim) = shiftlocus(Center(), convertmode(Projected, dim))
