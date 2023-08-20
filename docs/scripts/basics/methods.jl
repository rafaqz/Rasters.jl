# ## Methods that change the resolution or extent of an object

# Click through to the function documentation for more in-depth descriptions and examples.
#    
# | <div style="width:120px">Methods</div> |              Description                                        |
# | :------------------------ | :--------------------------------------------------------------------------- |
# | [`aggregate`](@ref)       | aggregate data by the same or different amounts for each axis.               |
# | [`disaggregate`](@ref)    | similarly disaggregate data.                                                 |
# | [`mosaic`](@ref)          | join rasters covering different extents into a single array or file.         |
# | [`crop`](@ref)            | shrink objects to specific dimension sizes or the extent of another object.  |
# | [`extend`](@ref)          | extend objects to specific dimension sizes or the extent of another object.  |
# | [`trim`](@ref)            | trims areas of missing values for arrays and across stack layers.            |
# | [`resample`](@ref)        | resample data to a different size and projection, or snap to another object. |
# | [`warp`](@ref)            | use `gdalwarp` on any object, e.g. a multidimensional NetCDF stack.          |
    
    
# ## Methods that change an objects values
#    
# !!! info
#       Note that most regular Julia methods, such as `replace`, work as for a standard
#       `Array`. These additional methods are commonly required in GIS applications.
#    
# | <div style="width:120px">Methods</div> |              Description                                        |
# | :------------------------ | :--------------------------------------------------------------------------- |
# | [`classify`](@ref)        | classify values into categories.                                             |
# | [`mask`](@ref)            | mask an object by a polygon or `Raster` along `X/Y`, or other dimensions.    |
# | [`replace_missing`](@ref) | replace all missing values in an object and update `missingval`.             |
#    
#    
# ## Point, polygon and table operation
#    
# | <div style="width:120px">Methods</div> |              Description                                        |
# | :------------------------ | :--------------------------------------------------------------------------- |
# | [`rasterize`](@ref)       | rasterize points and geometries.                                             |
# | [`extract`](@ref)         | extract values from points or geometries.                                    |
# | [`zonal`](@ref)           | calculate zonal statistics for an object masked by geometries.               |    
#    
# ## Methods to load, write and modify data sources
#    
# | <div style="width:120px">Methods</div> |              Description                                   |
# | :------------------------ | :---------------------------------------------------------------------- |
# | [`modify`](@ref)          | replace the data in objects. Useful to e.g. move objects to/from a GPU. |
# | [`read`](@ref)            | read data to memory if it is on disk.                                   |
# | [`read!`](@ref)           | read data to predefined memory.                                         |
# | [`open`](@ref)            | open the underlying data for manually reading or writing.               |
# | [`write`](@ref)           | write objects to file.                                                  |
    
