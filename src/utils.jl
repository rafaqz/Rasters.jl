filter_ext(path, ext) = filter(filename -> splitext(filename)[2] == ext, readdir(path))

front(s::AbstractString) = s[1:end-1];

windowsize(window::Tuple{Int,Vararg}) = (windowsize(tail(window))...,)
windowsize(window::Tuple) = (length(window[1]), windowsize(tail(window))...)
windowsize(window::Tuple{}) = ()

maybewindow2indices(A, window::Tuple) =
    maybewindow2indices(A, dims(A), window::Tuple)
maybewindow2indices(A, dims::Tuple, window::Tuple) =
    window == () ? () : to_indices(A, dims2indices(dims, window))

readwindowed(A::AbstractGeoArray, window::Tuple{}) = GeoArray(A)
readwindowed(A, window::Tuple{}) = Array(A)
readwindowed(A, window::Tuple{}, I...) = A[I...]
readwindowed(A, window::Tuple, I...) = A[Base.reindex(window, I)...]
readwindowed(A, window::Tuple) = A[window...]
readwindowed(A) = Array(A)
readwindowed(A, I...) = A[I...]

@inline getmeta(A::AbstractGeoArray, key, fallback) = getmeta(metadata(A), key, fallback)
@inline getmeta(m::Nothing, key, fallback) = fallback
@inline getmeta(m::Union{NamedTuple,Dict}, key, fallback) = key in keys(m) ?  m[key] : fallback
@inline getmeta(m::Metadata, key, fallback) = getmeta(val(m), key, fallback)

checkarrayorder(A, order::Order) = map(d -> checkarrayorder(d, order), dims(A))
checkarrayorder(A, order::Tuple) = map(checkarrayorder, dims(A), order)
checkarrayorder(dim::Dimension, order::Order) = begin
    ao = arrayorder(dim)
    ao == order || @warn "Array order for `$(DD.basetypeof(order))` is `$ao`, usualy `$order`"
end


cleankeys(keys) = Tuple(Symbol.(keys))
