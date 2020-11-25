abstract type GeoXDim{T,Mo,Me} <: XDim{T,Mo,Me} end
abstract type GeoYDim{T,Mo,Me} <: YDim{T,Mo,Me} end
abstract type GeoZDim{T,Mo,Me} <: ZDim{T,Mo,Me} end

"""
    Lon <: XDim <: Dimension
    Lon(val=:)

Longitude [`Dimension`]($DDdimdocs).

## Example:
```julia
longdim = Lon(10:10:100)
# Or
val = A[Lon(1)]
# Or
mean(A; dims=Lon)
```
"""
@dim Lon GeoXDim "Longitude"

"""
    Lat <: YDim <: Dimension
    Lat(val=:)

Latitude [`Dimension`]($DDdimdocs).

## Example:
```julia
vertdim = Lat(10:10:100)
# Or
val = A[Lat(1)]
# Or
mean(A; dims=Lat)
```
"""
@dim Lat GeoYDim "Latitude"

"""
    Vert <: ZDim <: Dimension
    Vert(val=:)

Vertical [`Dimension`]($DDdimdocs).

## Example:
```julia
vertdim = Vert(10:10:100)
# Or
val = A[Vert(1)]
# Or
mean(A; dims=Vert)
```
"""
@dim Vert GeoZDim "Vertical"

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

Get the bounds converted to the `usercrs` value.

Whithout ArchGDAL loaded, this is just the regular bounds.
"""
function mappedbounds end

mappedbounds(dims::Tuple) = map(mappedbounds, dims)
mappedbounds(dim::Dimension) = bounds(dim)
mappedbounds(dim::Union{Lat,Lon}) = mappedbounds(mode(dim), dim)
mappedbounds(::Mapped, dim) = bounds(dim)
@noinline mappedbounds(mode::IndexMode, dim) =
    if mode isa Projected
        error("Load ArchGDAL to convert Projected mode bounds to mapped")
    else
        error("cannot get mapped bounds of a $(nameof(typeof(mode))) mode dim")
    end


projectedbounds(dims::Tuple) = map(projectedbounds, dims)
projectedbounds(dim::Dimension) = bounds(dim)
projectedbounds(dim::Union{Lat,Lon}) = projectedbounds(mode(dim), dim)
projectedbounds(::Projected, dim) = bounds(dim)
@noinline projectedbounds(mode::IndexMode, dim) =
    if mode isa Mapped
        error("Load ArchGDAL to convert Mapped mode dim to projected")
    else
        error("cannot get projected bounds of a $(nameof(typeof(mode))) mode dim")
    end


"""
    mappedindex(x)

Get the index value of a dimension converted to the `usercrs` value.

Whithout ArchGDAL loaded, this is just the regular dim value.
"""
function mappedindex end

mappedindex(dims::Tuple) = map(mappedindex, dims)
mappedindex(dim::Dimension) = index(dim)
mappedindex(dim::Union{Lat,Lon}) = mappedindex(mode(dim), dim)
mappedindex(::Mapped, dim) = index(dim)
@noinline mappedindex(mode::IndexMode, dim) =
    if mode isa Projected
        error("Load ArchGDAL to convert Projected mode index to mapped")
    else
        error("cannot get mapped index of a $(nameof(typeof(mode))) mode dim")
    end

projectedindex(dims::Tuple) = map(projectedindex, dims)
projectedindex(dim::Dimension) = index(dim)
projectedindex(dim::Union{Lat,Lon}) = projectedindex(mode(dim), dim)
projectedindex(::Projected, dim) = index(dim)
@noinline projectedindex(mode::IndexMode, dim) =
    if mode isa Mapped
        error("Load ArchGDAL to convert Mapped mode index to projected")
    else
        error("cannot get projected index of a $(nameof(typeof(mode))) mode dim")
    end

