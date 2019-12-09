filter_ext(path, ext) = filter(filename -> splitext(filename)[2] == ext, readdir(path))

front(s::AbstractString) = s[1:end-1];


applywindow(s::AbstractGeoStack, key::Key, I) = applywindow(windoworempty(s), I)
applywindow(s::AbstractGeoArray, I) = applywindow(windoworempty(s), I)
applywindow(dims::Tuple{}, I::Tuple) = I
applywindow(dims::Tuple, I::Tuple) = (applywindow(dims[1], I[1]), applywindow(tail(dims), tail(I))...)
applywindow(dims::Tuple{}, I::Tuple{}) = ()  
applywindow(dims::Colon, I) = I  
applywindow(dims, I) = dims[I]  

window2indices(x) = window2indices(dims(x), window(x))
window2indices(x, window) = window 
window2indices(dims, window::Tuple{Vararg{<:AbstractDimension}}) = dims2indices(dims, window)

windoworempty(x) = (w = window(x); w == () ? () : window2indices(x))

windowsize(window::Tuple, sze::Tuple) = map(windowsize, window, sze)
windowsize(window::Colon, sze) = sze
windowsize(window, sze) = window
