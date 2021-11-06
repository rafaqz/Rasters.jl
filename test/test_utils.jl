# Loader for external sources
maybedownload(url, filename=splitdir(url)[2]) = begin
    dirpath = joinpath(dirname(pathof(Rasters)), "../test/data")
    mkpath(dirpath)
    filepath = joinpath(dirpath, filename) 
    isfile(filepath) || download(url, filepath)
    filepath
end
