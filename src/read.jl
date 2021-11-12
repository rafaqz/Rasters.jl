"""
    read(A::AbstractRaster)
    read(A::AbstractRasterStack)
    read(A::AbstractRasterSeries)

`read` will move a Rasters.jl object completely to memory.
"""
function Base.read(x::Union{AbstractRaster,AbstractRasterStack,AbstractRasterSeries})
    modify(x) do ds
        # Some backends don't implement `Array` properly.
        Array{eltype(ds),ndims(ds)}(undef, size(ds)) .= ds
    end
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
    src = Raster(filename;
        dims=dims(dst), refdims=refdims(dst), name=name(dst),
        metadata=metadata(dst), missingval=missingval(dst),
    )
    read!(src, dst)
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
    src = RasterStack(filenames;
        dims=dims(dst), refdims=refdims(dst), keys=keys(dst), metadata=metadata(dst),
        layermetadata=DD.layermetadata(dst), missingval=missingval(dst),
    )
    read!(src, dst)
end
