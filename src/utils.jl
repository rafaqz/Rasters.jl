filter_ext(path, ext::AbstractString) = filter(fn -> splitext(fn)[2] == ext, readdir(path))
filter_ext(path, exts::Union{Tuple,AbstractArray}) = 
    filter(fn -> splitext(fn)[2] in exts, readdir(path))
filter_ext(path, ext::Nothing) = readdir(path)

maybewindow2indices(A, window::Tuple) = maybewindow2indices(A, dims(A), window::Tuple)
maybewindow2indices(A, dims::Tuple, window::Tuple) =
    window == () ? () : to_indices(A, DD.dims2indices(dims, window))

# Read from the paraent dataset, using the window indices if they exist
# We try to load as little data from disk as possible.
readwindowed(A::AbstractGeoArray, window::Tuple{}) = GeoArray(A)
readwindowed(A, window::Tuple{}) = Array(A)
readwindowed(A, window::Tuple{}, I...) = A[I...]
readwindowed(A, window::Tuple, I...) = A[Base.reindex(window, I)...]
readwindowed(A, window::Tuple) = readwindowed(A, window...)
readwindowed(A, i, I...) = A[i, I...]
readwindowed(A) = Array(A)

# Get a metadata field
getmeta(A::AbstractGeoArray, key, fallback) = getmeta(metadata(A), key, fallback)
getmeta(m::Metadata, key, fallback) = get(val(m), key, fallback)
getmeta(m::NoMetadata, key, fallback) = fallback

# Check that array order matches expectation
checkarrayorder(A, order::Order) = map(d -> checkarrayorder(d, order), dims(A))
checkarrayorder(A, order::Tuple) = map(checkarrayorder, dims(A), order)
checkarrayorder(dim::Dimension, order::Order) =
    arrayorder(dim) == order || @warn "Array order for `$(DD.basetypeof(order))` is `$(arrayorder(dim))`, usualy `$order`"

checkindexorder(A, order::Order) = map(d -> checkindexorder(d, order), dims(A))
checkindexorder(A, order::Tuple) = map(checkindexorder, dims(A), order)
checkindexorder(dim::Dimension, order::Order) =
    indexorder(dim) == order || @warn "Array order for `$(DD.basetypeof(order))` is `$(indexorder(dim))`, usualy `$order`"

cleankeys(keys) = Tuple(Symbol.(keys))

function _to_sampled(A) 
    if any(m -> !(m isa AbstractSampled), mode(A))
        rebuild(A; data=data(A), dims=map(_to_sampled, dims(A)))
    end
end
_to_sampled(mode, dim) = dim
_to_sampled(mode::AbstractSampled, dim::SpatialDim) = dim
_to_sampled(mode, dim::SpatialDim) = rebuild(dim; mode=Sampled(mode(dim)))
