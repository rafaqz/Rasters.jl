"""
    read(A::AbstractGeoArray)
    read(A::AbstractGeoStack)
    read(A::AbstractGeoSeries)

`read` will move a GeoData.jl object completely to memory, or make
a copy if already in memory.
```
"""
Base.read(x::Union{AbstractGeoArray,AbstractGeoStack,AbstractGeoSeries}) = modify(Array, x)

"""
    read!(filename, A::AbstractGeoArray)
    read!(filename, A::AbstractGeoStack)
    read!(filename, A::AbstractGeoSeries)

`read!` will move a GeoData.jl object completely to memory, or make
a copy if already in memory.
```
"""
Base.read!(src::AbstractGeoArray, dst::AbstractArray) = dst .= src
function Base.read!(
    src::AbstractGeoStack, dst::AbstractGeoStack{Union{NamedTuple{Keys},FileStack{<:Any,Keys}}}
) where Keys
    map(Keys) do k
        read!(dst[k], src[k])
    end
    return dst
end
function Base.read!(src::AbstractGeoSeries, dst::AbstractGeoSeries)
    map(read!, src, dst)
    return dst
end

function Base.read!(filename::AbstractString, dst::AbstractGeoArray)
    src = geoarray(filename;
        dims=dims(dst),
        refdims=refdims(dst),
        metadata=metadata(dst),
        missingval=missingval(dst),
    )
    read!(src, dst)
end
function Base.read!(filename::AbstractString, dst::AbstractGeoStack)
    src = stack(filename;
        dims=dims(dst),
        refdims=refdims(dst),
        metadata=metadata(dst),
        missingval=missingval(dst),
    )
    read!(src, dst)
end
function Base.read!(filenames::AbstractVector{<:AbstractString}, dst::AbstractGeoStack)
    src = stack(filename;
        dims=dims(dst),
        refdims=refdims(dst),
        metadata=metadata(dst),
        missingval=missingval(dst),
    )
    read!(src, dst)
end
function Base.read!(filenames::AbstractVector{<:AbstractString}, dst::AbstractGeoSeries)
    map((fn, d) -> read!(fn, d), filenames, dst)
    return dst
end
