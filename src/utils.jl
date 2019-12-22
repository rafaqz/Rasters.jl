filter_ext(path, ext) = filter(filename -> splitext(filename)[2] == ext, readdir(path))

front(s::AbstractString) = s[1:end-1];

windowsize(window::Tuple{Int,Vararg}) = (windowsize(tail(window))...,)
windowsize(window::Tuple) = (length(window[1]), windowsize(tail(window))...)
windowsize(window::Tuple{}) = ()

maybewindow2indices(A, dims, window) = 
    window == () ? () : to_indices(A, dims2indices(dims, window))

readwindowed(A, window::Tuple{}) = (println(typeof(A)); Array(A))
readwindowed(A, window::Tuple{}, I...) = A[I...]
readwindowed(A, window, I...) = A[Base.reindex(window, I)...]
readwindowed(A, window) = A[window...]
