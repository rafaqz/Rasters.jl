# Loader for external sources
geturl(url, filename=splitdir(url)[2]) = begin
    dirpath = joinpath(dirname(pathof(GeoData)), "../test/data")
    mkpath(dirpath)
    filepath = joinpath(dirpath, filename) 
    isfile(filepath) || download(url, filepath)
    filepath
end
