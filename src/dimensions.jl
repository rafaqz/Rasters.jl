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
@dim Lon XDim "Longitude" "Lon"

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
@dim Lat YDim "Latitude" "Lat"

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
@dim Vert ZDim "Vertical" "Vert"

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
