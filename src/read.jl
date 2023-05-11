"""
    read(A::AbstractRaster)
    read(A::AbstractRasterStack)
    read(A::AbstractRasterSeries)

`read` will move a Rasters.jl object completely to memory.
"""
function Base.read(x::Union{AbstractRaster,AbstractRasterStack,AbstractRasterSeries})
    modify(Array, x)
end

"""
    read!(src::Union{AbstractString,AbstractRaster}, dst::AbstractRaster)
    read!(src::Union{AbstractString,AbstractRasterStack}, dst::AbstractRasterStack)
    read!(scr::AbstractRasterSeries, dst::AbstractRasterSeries)

`read!` will copy the data from `src` to the object `dst`. 

`src` can be an object or a file-path `String`.
"""
Base.read!(src::AbstractRaster, dst::AbstractArray) = dst .= src
function Base.read!(src::AbstractRasterStack, dst::AbstractRasterStack)
    map(k -> read!(src[k], dst[k]), keys(dst))
    return dst
end
function Base.read!(src::AbstractRasterSeries, dst::AbstractRasterSeries)
    map(read!, src, dst)
    return dst
end

# Filename methods
function Base.read!(filename::AbstractString, dst::AbstractRaster)
    read!(Raster(filename; lazy=true), dst)
end
function Base.read!(filenames::Union{NamedTuple,<:AbstractVector{<:AbstractString}}, dst::AbstractRasterStack)
    _readstack!(filenames, dst)
end
function Base.read!(filenames::AbstractString, dst::AbstractRasterStack)
    _readstack!(filenames, dst)
end
function Base.read!(filenames::AbstractVector{<:Union{AbstractString,NamedTuple}}, dst::AbstractRasterSeries)
    map((fn, d) -> read!(fn, d), filenames, dst)
    return dst
end

function _readstack!(filenames, dst)
    read!(RasterStack(filenames; lazy=true), dst)
end
