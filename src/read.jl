"""
    read(A::AbstractRaster)
    read(A::AbstractRasterStack)
    read(A::AbstractRasterSeries)

`read` will move a Rasters.jl object completely to memory.
"""
function Base.read(x::Union{AbstractRaster,AbstractRasterStack,AbstractRasterSeries})
    _checkmem(x)
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

_sizeof(A::AbstractArray{T}) where T = sizeof(T) * prod(size(A))
_sizeof(st::AbstractRasterStack) = sum(_sizeof, layers(st))
_sizeof(s::AbstractRasterSeries) =
    length(s) == 0 ? 0 : _sizeof(first(s)) * prod(size(s))

function _checkmem(x)
    required_mem = _sizeof(x)
    Sys.free_memory() > required_mem || _no_mem_error(required_mem)
end

_no_mem_error(required_mem) =
    error("required memory $required_mem is greater than system memory $(Sys.free_memory()). Use `lazy=true` if you are loading dataset, and only call `read` on a subset after `view`.")
