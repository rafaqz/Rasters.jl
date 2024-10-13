#=
# `resample` tutorial - warping a raster
=#
using Rasters, ArchGDAL

#=
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

