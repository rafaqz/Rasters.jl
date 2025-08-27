
const SpatialDim = Union{XDim,YDim,ZDim}

"""
    Band <: Dimension

    Band(val=:)

Band [`Dimension`]($DDdimdocs) for multi-band rasters.

## Example
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
    Lat <: Dimension

    Lat(val=:)

Used for holding degrees north lookups.

Will error on lookup construction if metadata of `units="degrees_north"` does not exist.
"""
@dim Lat

"""
    Lon <: Dimension

    Lon(val=:)

Used for holding degrees east lookups.

Will error on lookup construction if metadata of `units="degrees_east"` does not exist.
"""
@dim Lon

"""
    Geometry <: Dimension

    Geometry(geoms)

Geometry [`Dimension`]($DDdimdocs) for vector data cubes.

## Example
```julia
geomdim = Geometry(GeometryLookup(polygons))
# Or
val = A[Geometry(1)]
# Or
val = A[Geometry(Touches(other_geom))] # this is automatically accelerated by spatial trees!
# Or
mean(A; dims=Geometry)
```
"""
@dim Geometry
