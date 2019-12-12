# Loader for external sources
geturl(url, filename=splitdir(url)[2]) = begin
    isfile(filename) || download(url, filename)
    filename
end
