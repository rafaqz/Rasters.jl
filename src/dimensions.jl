
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
    Geometry <: Dimension
    Geometry(geoms)

Geometry [`Dimension`]($DDdimdocs) for vector data cubes.

"""
@dim Geometry


# Below is the example docstring for the Geometry dimension, which should be added back to the docs
# when geometry lookups are fully implemented.
#=
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
=#