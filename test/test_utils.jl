# Loader for external sources
geturl(url, filename=splitdir(url)[2]) = begin
    filepath = joinpath(dirname(pathof(GeoData)), "../test/data", filename) 
    println(filepath)
    isfile(filepath) || download(url, filepath)
    filepath
end
