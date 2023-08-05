# Most base methods work as for regular julia `Array`s, such as `reverse` and
# rotations like `rotl90`. Base, statistics and linear algebra methods like `mean`
# that take a `dims` argument can also use the dimension name. 

# ## Mean over the time dimension:

using Rasters, Statistics, RasterDataSources

A = Raster(WorldClim{BioClim}, 5)

# Then we do the meand over the `X` dimension

mean(A, dims=X) # Ti if time were available would also be possible

# `broadcast` works lazily from disk when `lazy=true`, and is only applied when data
# is directly indexed. Adding a dot to any function will use broadcast over a `Raster`
# just like an `Array`. 

# ## Broadcasting

# For a disk-based array `A`, this will only be applied when indexing occurs or
# when we [`read`](@ref) the array.
 
A .*= 2

# To broadcast directly to disk, we need to open the file in write mode:

# ```
# open(Raster(filename); write=true) do O
#    O .*= 2
# end
# ```

# To broadcast over a `RasterStack` use `map`, which applies a function to
# the raster layers of the stack.

# ```
# newstack = map(stack) do raster
#    raster .* 2
# end
# ```

# ## Modifying object properties

# `rebuild` can be used to modify the fields of an object, generating a new object
# (but possibly holding the same arrays or files).

# If you know that a file had an incorrectly specified missing value, you can do:

# ### rebuild

A = Raster(WorldClim{BioClim}, 5)
rebuild(A; missingval=-9999)

# (`replace_missing` will actualy replace the current values).

# Or if you need to change the name of the layer:

# Then use `rebuild` as
B = rebuild(A; name=:temperature)
B.name

# ### replace_missing

replace_missing(A, missingval=-9999)

# ### set

# `set` can be used to modify the nested properties of an objects dimensions, that
# are more difficult to change with `rebuild`. `set` works on the principle that
# dimension properties can only be in one specific field, so we generally don't
# have to specify which one it is. `set` will also try to update anything affected
# by a change you make.

# This will set the `X` axis to specify points, instead of intervals:

using Rasters: Points
set(A, X => Points)

# We can also reassign dimensions, here `X` becomes `Z`:

set(A, X => Z)

# `setcrs(A, crs)` and `setmappedcrs(A, crs)` will set the crs value/s of an
# object to any `GeoFormat` from GeoFormatTypes.jl.