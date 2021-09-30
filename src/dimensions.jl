
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
