
const SpatialDim = Union{XDim,YDim,ZDim}

"""
    Band <: Dimension

    Band(val=:)

Band [`Dimension`]($DDdimdocs) for multi-band rasters.

## Example:
```julia
banddim = Band(10:10:100)
# Or
val = A[Band(1)]
# Or
mean(A; dims=Band)
```
"""
@dim Band

"""
    mappedbounds(x)

Get the bounds converted to the [`mappedcrs`](@ref) value.

Whithout ArchGDAL loaded, this is just the regular bounds.
"""
function mappedbounds end

mappedbounds(dims::Tuple) = map(mappedbounds, dims)
mappedbounds(dim::Dimension) = mappedbounds(mode(dim), dim)
mappedbounds(::IndexMode, dim) = bounds(dim)
mappedbounds(mode::Projected, dim) = mappedbounds(mappedcrs(mode), mode, dim)
mappedbounds(mappedcrs::Nothing, mode::Projected, dim) =
    error("No mappedcrs attached to $(name(dim)) dimension")
mappedbounds(mappedcrs::GeoFormat, mode::Projected, dim) =
    _sort(reproject(crs(mode), mappedcrs, dim, bounds(dim)))

projectedbounds(dims::Tuple) = map(projectedbounds, dims)
projectedbounds(dim::Dimension) = projectedbounds(mode(dim), dim)
projectedbounds(::IndexMode, dim) = bounds(dim)
projectedbounds(mode::Mapped, dim) = projectedbounds(crs(mode), mode, dim)
projectedbounds(crs::Nothing, mode::Mapped, dim) =
    error("No projection crs attached to $(name(dim)) dimension")
projectedbounds(crs::GeoFormat, mode::Mapped, dim) =
    _sort(reproject(mappedcrs(mode), crs, dim, bounds(dim)))

_sort((a, b)) = a <= b ? (a, b) : (b, a)

"""
    mappedindex(x)

Get the index value of a dimension converted to the `mappedcrs` value.

Whithout ArchGDAL loaded, this is just the regular dim value.
"""
function mappedindex end

mappedindex(dims::Tuple) = map(mappedindex, dims)
mappedindex(dim::Dimension) = _mappedindex(mode(dim), dim)

_mappedindex(::IndexMode, dim::Dimension) = index(dim)
_mappedindex(mode::Projected, dim::Dimension) = _mappedindex(mappedcrs(mode), mode, dim)
_mappedindex(mappedcrs::Nothing, mode::Projected, dim) =
    error("No mappedcrs attached to $(name(dim)) dimension")
_mappedindex(mappedcrs::GeoFormat, mode::Projected, dim) =
    reproject(crs(dim), mappedcrs, dim, index(dim))

projectedindex(dims::Tuple) = map(projectedindex, dims)
projectedindex(dim::Dimension) = _projectedindex(mode(dim), dim)

_projectedindex(::IndexMode, dim::Dimension) = index(dim)
_projectedindex(mode::Mapped, dim::Dimension) = _projectedindex(crs(mode), mode, dim)
_projectedindex(crs::Nothing, mode::Mapped, dim::Dimension) =
    error("No projection crs attached to $(name(dim)) dimension")
_projectedindex(crs::GeoFormat, mode::Mapped, dim::Dimension) =
    reproject(mappedcrs(dim), crs, dim, index(dim))
