using ArchGDAL

# Loader for external sources
maybedownload(url, filename=splitdir(url)[2]) = begin
    dirpath = joinpath(dirname(pathof(Rasters)), "../test/data")
    mkpath(dirpath)
    filepath = joinpath(dirpath, filename) 
    isfile(filepath) || download(url, filepath)
    filepath
end

# create `N` random rasters with eltype `type` and x,y,band size `size`
function temporary_random_rasters(f, N, size, type=UInt8)
    filenames = [tempname() * ".tif" for _ in 1:N]
    try
        for f in filenames
            write(f, Raster(rand(type, size); dims=(X(1:size[1]), Y(1:size[2]), Band(1:size[3]))))
        end
        f(filenames)
    finally
        rm.(filenames; force=true)
    end
end