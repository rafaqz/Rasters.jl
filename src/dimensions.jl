
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

[`Dimension`]($DDdimdocs) for a latitude axis (`degrees_north` in CF terms).

CF (Climate and Forecast) datasets with 2D coordinate variables expose
their latitude as a `Lat`-typed inner dim of a `ProjectedArrayLookup`. The
CF role - latitude vs longitude - is decided at load time by the
`_unaligned_lookup` reader in `sources/commondatamodel.jl` based on the
`units` / `standard_name` attributes of the coordinate variable.
"""
@dim Lat

"""
    Lon <: Dimension

    Lon(val=:)

[`Dimension`]($DDdimdocs) for a longitude axis (`degrees_east` in CF terms).

See [`Lat`](@ref) for how the CF role is decided at load time.
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
