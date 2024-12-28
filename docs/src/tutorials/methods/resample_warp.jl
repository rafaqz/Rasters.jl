#=
# `resample` tutorial - warping a raster

```@meta
CollapsedDocStrings=true
```
=#
using Rasters, ArchGDAL
using RasterDataSources
using NaNStatistics
using CairoMakie
#=
## What is resampling?

**[`resample`](@ref)** "re-samples" the 
data by interpolation and can also aggregate or disaggregate, changing the resolution.
It always returns a `Regular` lookup (like a range), and is the most flexible of the 
resampling methods.

This uses GDAL's `gdalwarp` algorithm under the hood.  You can call that via [`warp`](@ref)
if you need more control, but generally `resample` is sufficient. 

Rasters.jl has a few other methods to change the lookups of a raster.  These are:
- [`reproject`](@ref), which simply directly reprojects the lookup axes 
  (but is **only usable for specific cases**, where the source and destination 
  coordinate systems are both cylindrical, like the long-lat, Mercator, or Web-Mercator projections.) 

  This is a lossless operation and keeps the data exactly the same - only the axes are changed. 

- [`aggregate`](@ref) and [`disaggregate`](@ref), which change the resolution of 
  the raster by clumping ([`aggregate`](@ref)) or splitting ([`disaggregate`](@ref)) cells.

  They can't change cells fractionally, and can't change the projection or coordinate system.

Of all these methods, **`resample`** is the most flexible and powerful, and is what we will focus on here.  
It is, however, also the slowest.  So if another method is applicable to your problem, you should consider it.

## How `resample` works

`resample` uses GDAL's `gdalwarp` algorithm under the hood.  This is a battle-tested algorithm
and is generally pretty robust.  However, it has the following limitations:
- It always assumes cell-based sampling, instead of point-based sampling.  This does mean that 
  point-based rasters are converted to cell-based sampling.
- It can only accept some primitive types for the input data, since that data is passed directly to a C library.
  Things like `RGB` or user-defined types are not usually supported.

`resample` allows you to specify several methods, which you can see if you expand the docstring below.

```@docs; canonical=false
resample
```







Topics:
- What is resampling?
    - When to resample vs reproject
    - Things to keep in mind
        - GDAL always changes the locus to cell sampling, you can reset this by using `shiftlocus`
        - You can in fact resample to another raster, if you want perfect alignment.
            - This doesn't work for irregularly sampled rasters.
        - Resampling is a lossy operation and takes time.  Try to avoid repeatedly resampling, and if you must do so, crop or trim the raster as much as you can first.
- Show the different resampling methods, in a grid
- Show some different projections and ways of constructing them
- Show how to use `size` and `res` to change the resolution of a raster
- Show how to use `warp` to reproject a raster
=#

