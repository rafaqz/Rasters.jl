# Loader for external sources
geturl(url, filename=splitdir(url)[2]) = begin
    filepath = joinpath(dirname(pathof(GeoData)), "../test", filename) 
    println(filepath)
    isfile(filepath) || download(url, filepath)
    filepath
end
