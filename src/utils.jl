filter_ext(path, ext) = filter(filename -> splitext(filename)[2] == ext, readdir(path))

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

getmeta(A::AbstractGeoArray, key, fallback) = getmeta(metadata(A), key, fallback)
getmeta(m::Nothing, key, fallback) = fallback
getmeta(m::Union{NamedTuple,Dict}, key, fallback) = key in keys(m) ?  m[key] : fallback
getmeta(m::Metadata, key, fallback) = getmeta(val(m), key, fallback)

checkarrayorder(A, order::Order) = map(d -> checkarrayorder(d, order), dims(A))
checkarrayorder(A, order::Tuple) = map(checkarrayorder, dims(A), order)
checkarrayorder(dim::Dimension, order::Order) = begin
    ao = arrayorder(dim)
    ao == order || @warn "Array order for `$(DD.basetypeof(order))` is `$ao`, usualy `$order`"
end

cleankeys(keys) = Tuple(Symbol.(keys))
