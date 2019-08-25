filter_ext(path, ext) = filter(filename -> splitext(filename)[2] == ext, readdir(path))

front(s::AbstractString) = s[1:end-1];
