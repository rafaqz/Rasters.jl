
## Point, polygon and table operations {#Point,-polygon-and-table-operations}

| Methods                                                 | Description                                                    |
|:------------------------------------------------------- |:-------------------------------------------------------------- |
| [`rasterize`](/api#Rasters.rasterize)                   | rasterize points and geometries.                               |
| [`extract`](/api#Rasters.extract)                       | extract values from points or geometries.                      |
| [`zonal`](/api#Rasters.zonal-Tuple{Any,%20RasterStack}) | calculate zonal statistics for an object masked by geometries. |


## Methods that change the resolution or extent of an object {#Methods-that-change-the-resolution-or-extent-of-an-object}

Click through to the function documentation for more in-depth descriptions and examples.

| Methods                                                                                                                                                 | Description                                                                  |
|:------------------------------------------------------------------------------------------------------------------------------------------------------- |:---------------------------------------------------------------------------- |
| [`aggregate`](/api#Rasters.aggregate)                                                                                                                   | aggregate data by the same or different amounts for each axis.               |
| [`disaggregate`](/api#Rasters.disaggregate)                                                                                                             | similarly disaggregate data.                                                 |
| [`mosaic`](/api#Rasters.mosaic-Tuple{Function,%20Union{AbstractRaster,%20AbstractRasterStack},%20Vararg{Union{AbstractRaster,%20AbstractRasterStack}}}) | join rasters covering different extents into a single array or file.         |
| [`crop`](/api#Rasters.crop)                                                                                                                             | shrink objects to specific dimension sizes or the extent of another object.  |
| [`extend`](/api#Rasters.extend)                                                                                                                         | extend objects to specific dimension sizes or the extent of another object.  |
| [`trim`](/api#Rasters.trim-Tuple{Union{AbstractRaster,%20AbstractRasterStack}})                                                                         | trims areas of missing values for arrays and across stack layers.            |
| [`resample`](/api#Rasters.resample-Tuple)                                                                                                               | resample data to a different size and projection, or snap to another object. |
| [`warp`](/api#Rasters.warp-Tuple)                                                                                                                       | use `gdalwarp` on any object, e.g. a multidimensional NetCDF stack.          |


## Methods that change an objects values {#Methods-that-change-an-objects-values}

| Methods                                                          | Description                                                               |
|:---------------------------------------------------------------- |:------------------------------------------------------------------------- |
| [`classify`](/api#Rasters.classify)                              | classify values into categories.                                          |
| [`mask`](/api#Rasters.mask-Tuple{Any})                           | mask an object by a polygon or `Raster` along `X/Y`, or other dimensions. |
| [`replace_missing`](/tutorials/array_operations#replace_missing) | replace all missing values in an object and update `missingval`.          |


::: tip Info

Note that most regular Julia methods, such as `replace`, work as for a standard   `Array`. These additional methods are commonly required in GIS applications.

:::

## Methods to load, write and modify data sources {#Methods-to-load,-write-and-modify-data-sources}

| Methods                                                                                              | Description                                                             |
|:---------------------------------------------------------------------------------------------------- |:----------------------------------------------------------------------- |
| [`modify`](/api#DimensionalData.modify-Tuple{Any,%20AbstractRasterSeries})                           | replace the data in objects. Useful to e.g. move objects to/from a GPU. |
| [`read`](/api#Base.read-Tuple{Union{AbstractRaster,%20AbstractRasterSeries,%20AbstractRasterStack}}) | read data to memory if it is on disk.                                   |
| [`read!`](/api#Base.read!-Tuple{AbstractRaster,%20AbstractArray})                                    | read data to predefined memory.                                         |
| [`open`](/api#Base.open-Tuple{Function,%20AbstractRaster})                                           | open the underlying data for manually reading or writing.               |
| [`write`](/api#Base.write-Tuple{AbstractString,%20AbstractRasterSeries})                             | write objects to file.                                                  |

