

create(filename, A::AbstractRaster{T}; kw...) where T = create(filename, T, A; kw...)
function create(filename, T, A::AbstractRaster; 
    name=name(A), metadata=metadata(A), missingval=missingval(A), kw...
)
    create(filename, T, dims(A); parent=parent(A), name, metadata, missingval, kw...)
end
function create(filename::AbstractString, T::Type, dims::Tuple; parent=nothing, suffix=nothing, kw...)
    filename = _maybe_add_suffix(filename, suffix)
    create(filename, _sourcetype(filename), T, dims; kw...)
end
function create(filename::Nothing, T::Type, dims::Tuple; parent=nothing, suffix=nothing, missingval, kw...)
    T = isnothing(missingval) ? T : promote_type(T, typeof(missingval))
    data = isnothing(parent) ? Array{T}(undef, dims) : similar(parent, T, size(dims))
    Raster(data, dims; missingval, kw...)
end

_maybe_add_suffix(filename::Nothing, suffix) = nothing
_maybe_add_suffix(filename::Nothing, suffix::Nothing) = nothing
_maybe_add_suffix(filename, suffix::Nothing) = filename
function _maybe_add_suffix(filename, suffix)
    base, ext = splitext(filename)
    if string(suffix) == ""
        filename
    else
        return string(base, "_", suffix, ext)
    end
end
