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

const Lon = X
const Lat = Y


"""
    mappedbounds(x)

Get the bounds converted to the [`mappedcrs`](@ref) value.

Whithout ArchGDAL loaded, this is just the regular bounds.
"""
function mappedbounds end

mappedbounds(dims::Tuple) = map(mappedbounds, dims)
mappedbounds(dim::Dimension) = bounds(dim)
mappedbounds(dim::Union{Y,X}) = mappedbounds(mode(dim), dim)
mappedbounds(::Mapped, dim) = bounds(dim)
@noinline mappedbounds(mode::IndexMode, dim) =
    if mode isa Projected
        error("Load ArchGDAL to convert Projected mode bounds to mapped")
    else
        error("cannot get mapped bounds of a $(nameof(typeof(mode))) mode dim")
    end


projectedbounds(dims::Tuple) = map(projectedbounds, dims)
projectedbounds(dim::Dimension) = bounds(dim)
projectedbounds(dim::Union{Y,X}) = projectedbounds(mode(dim), dim)
projectedbounds(::Projected, dim) = bounds(dim)
@noinline projectedbounds(mode::IndexMode, dim) =
    if mode isa Mapped
        error("Load ArchGDAL to convert Mapped mode dim to projected")
    else
        error("cannot get projected bounds of a $(nameof(typeof(mode))) mode dim")
    end


"""
    mappedindex(x)

Get the index value of a dimension converted to the `mappedcrs` value.

Whithout ArchGDAL loaded, this is just the regular dim value.
"""
function mappedindex end

mappedindex(dims::Tuple) = map(mappedindex, dims)
mappedindex(dim::Dimension) = index(dim)
mappedindex(dim::Union{Y,X}) = mappedindex(mode(dim), dim)
mappedindex(::Mapped, dim) = index(dim)
@noinline mappedindex(mode::IndexMode, dim) =
    if mode isa Projected
        error("Load ArchGDAL to convert Projected mode index to mapped")
    else
        error("cannot get mapped index of a $(nameof(typeof(mode))) mode dim")
    end

projectedindex(dims::Tuple) = map(projectedindex, dims)
projectedindex(dim::Dimension) = index(dim)
projectedindex(dim::Union{Y,X}) = projectedindex(mode(dim), dim)
projectedindex(::Projected, dim) = index(dim)
@noinline projectedindex(mode::IndexMode, dim) =
    if mode isa Mapped
        error("Load ArchGDAL to convert Mapped mode index to projected")
    else
        error("cannot get projected index of a $(nameof(typeof(mode))) mode dim")
    end
