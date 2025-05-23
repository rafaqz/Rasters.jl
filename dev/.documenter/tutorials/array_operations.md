
# Array operations {#Array-operations}

Most base methods work as for regular julia `Array`s, such as `reverse` and rotations like `rotl90`. Base, statistics and linear algebra methods like `mean` that take a `dims` argument can also use the dimension name. 

## Mean over the time dimension: {#Mean-over-the-time-dimension:}

```julia
using Rasters, Statistics, RasterDataSources

A = Raster(WorldClim{BioClim}, 5)
```


```
┌ 2160×1080 Raster{Union{Missing, Float32}, 2} bio5 ┐
├───────────────────────────────────────────────────┴──────────────────── dims ┐
  ↓ X Projected{Float64} -180.0:0.16666666666666666:179.83333333333331 ForwardOrdered Regular Intervals{Start},
  → Y Projected{Float64} 89.83333333333333:-0.16666666666666666:-90.0 ReverseOrdered Regular Intervals{Start}
├──────────────────────────────────────────────────────────────────── metadata ┤
  Metadata{Rasters.GDALsource} of Dict{String, Any} with 1 entry:
  "filepath" => "./WorldClim/BioClim/wc2.1_10m_bio_5.tif"
├────────────────────────────────────────────────────────────────────── raster ┤
  missingval: missing
  extent: Extent(X = (-180.0, 179.99999999999997), Y = (-90.0, 90.0))
  crs: GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.25722...
└──────────────────────────────────────────────────────────────────────────────┘
    ↓ →    89.8333    89.6667    89.5       …  -89.6667  -89.8333  -90.0
 -180.0      missing    missing    missing     -15.399   -13.805   -14.046
 -179.833    missing    missing    missing     -15.9605  -14.607   -14.5545
    ⋮                                       ⋱                        ⋮
  179.5      missing    missing    missing     -18.2955  -16.7583  -16.72
  179.667    missing    missing    missing     -18.2847  -16.7513  -16.72
  179.833    missing    missing    missing  …  -17.1478  -15.4243  -15.701
```


Then we do the mean over the `X` dimension

```julia
mean(A, dims=X) # Ti if time were available would also be possible
```


```
┌ 1×1080 Raster{Union{Missing, Float32}, 2} bio5 ┐
├────────────────────────────────────────────────┴─────────────────────── dims ┐
  ↓ X Projected{Float64} -180.0:360.0:-180.0 ForwardOrdered Regular Intervals{Start},
  → Y Projected{Float64} 89.83333333333333:-0.16666666666666666:-90.0 ReverseOrdered Regular Intervals{Start}
├──────────────────────────────────────────────────────────────────── metadata ┤
  Metadata{Rasters.GDALsource} of Dict{String, Any} with 1 entry:
  "filepath" => "./WorldClim/BioClim/wc2.1_10m_bio_5.tif"
├────────────────────────────────────────────────────────────────────── raster ┤
  missingval: missing
  extent: Extent(X = (-180.0, 180.0), Y = (-90.0, 90.0))
  crs: GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.25722...
└──────────────────────────────────────────────────────────────────────────────┘
    ↓ →  89.8333    89.6667    89.5       …  -89.6667  -89.8333  -90.0
 -180.0    missing    missing    missing     -23.0768  -22.9373  -22.0094
```


`broadcast` works lazily from disk when `lazy=true`, and is only applied when data is directly indexed. Adding a dot to any function will use broadcast over a `Raster` just like an `Array`. 

## Broadcasting {#Broadcasting}

For a disk-based array `A`, this will only be applied when indexing occurs or when we [`read`](/api#Base.read-Tuple{Union{AbstractRaster,%20AbstractRasterSeries,%20AbstractRasterStack}}) the array.

```julia
A .*= 2
```


```
┌ 2160×1080 Raster{Union{Missing, Float32}, 2} bio5 ┐
├───────────────────────────────────────────────────┴──────────────────── dims ┐
  ↓ X Projected{Float64} -180.0:0.16666666666666666:179.83333333333331 ForwardOrdered Regular Intervals{Start},
  → Y Projected{Float64} 89.83333333333333:-0.16666666666666666:-90.0 ReverseOrdered Regular Intervals{Start}
├──────────────────────────────────────────────────────────────────── metadata ┤
  Metadata{Rasters.GDALsource} of Dict{String, Any} with 1 entry:
  "filepath" => "./WorldClim/BioClim/wc2.1_10m_bio_5.tif"
├────────────────────────────────────────────────────────────────────── raster ┤
  missingval: missing
  extent: Extent(X = (-180.0, 179.99999999999997), Y = (-90.0, 90.0))
  crs: GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.25722...
└──────────────────────────────────────────────────────────────────────────────┘
    ↓ →    89.8333    89.6667    89.5       …  -89.6667  -89.8333  -90.0
 -180.0      missing    missing    missing     -30.798   -27.61    -28.092
 -179.833    missing    missing    missing     -31.921   -29.214   -29.109
    ⋮                                       ⋱                        ⋮
  179.5      missing    missing    missing     -36.591   -33.5165  -33.44
  179.667    missing    missing    missing     -36.5695  -33.5025  -33.44
  179.833    missing    missing    missing  …  -34.2955  -30.8485  -31.402
```


To broadcast directly to disk, we need to open the file in write mode:

```julia
open(Raster(filename); write=true) do O
   O .*= 2
end
```


To broadcast over a `RasterStack` use `map`, which applies a function to the raster layers of the stack.

```julia
newstack = map(stack) do raster
   raster .* 2
end
```


## Modifying object properties {#Modifying-object-properties}

`rebuild` can be used to modify the fields of an object, generating a new object (but possibly holding the same arrays or files).

If you know that a file had an incorrectly specified missing value, you can do:

### rebuild {#rebuild}

```julia
A = Raster(WorldClim{BioClim}, 5)
rebuild(A; missingval=-9999)
```


```
┌ 2160×1080 Raster{Union{Missing, Float32}, 2} bio5 ┐
├───────────────────────────────────────────────────┴──────────────────── dims ┐
  ↓ X Projected{Float64} -180.0:0.16666666666666666:179.83333333333331 ForwardOrdered Regular Intervals{Start},
  → Y Projected{Float64} 89.83333333333333:-0.16666666666666666:-90.0 ReverseOrdered Regular Intervals{Start}
├──────────────────────────────────────────────────────────────────── metadata ┤
  Metadata{Rasters.GDALsource} of Dict{String, Any} with 1 entry:
  "filepath" => "./WorldClim/BioClim/wc2.1_10m_bio_5.tif"
├────────────────────────────────────────────────────────────────────── raster ┤
  missingval: -9999.0f0
  extent: Extent(X = (-180.0, 179.99999999999997), Y = (-90.0, 90.0))
  crs: GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.25722...
└──────────────────────────────────────────────────────────────────────────────┘
    ↓ →    89.8333    89.6667    89.5       …  -89.6667  -89.8333  -90.0
 -180.0      missing    missing    missing     -15.399   -13.805   -14.046
 -179.833    missing    missing    missing     -15.9605  -14.607   -14.5545
    ⋮                                       ⋱                        ⋮
  179.5      missing    missing    missing     -18.2955  -16.7583  -16.72
  179.667    missing    missing    missing     -18.2847  -16.7513  -16.72
  179.833    missing    missing    missing  …  -17.1478  -15.4243  -15.701
```


(`replace_missing` will actually replace the current values).

Or if you need to change the name of the layer:

Then use `rebuild` as

```julia
B = rebuild(A; name=:temperature)
B.name
```


```
:temperature
```


### replace_missing {#replace_missing}

```julia
replace_missing(A, missingval=-9999)
```


```
┌ 2160×1080 Raster{Float32, 2} bio5 ┐
├───────────────────────────────────┴──────────────────────────────────── dims ┐
  ↓ X Projected{Float64} -180.0:0.16666666666666666:179.83333333333331 ForwardOrdered Regular Intervals{Start},
  → Y Projected{Float64} 89.83333333333333:-0.16666666666666666:-90.0 ReverseOrdered Regular Intervals{Start}
├──────────────────────────────────────────────────────────────────── metadata ┤
  Metadata{Rasters.GDALsource} of Dict{String, Any} with 1 entry:
  "filepath" => "./WorldClim/BioClim/wc2.1_10m_bio_5.tif"
├────────────────────────────────────────────────────────────────────── raster ┤
  missingval: -9999.0f0
  extent: Extent(X = (-180.0, 179.99999999999997), Y = (-90.0, 90.0))
  crs: GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.25722...
└──────────────────────────────────────────────────────────────────────────────┘
    ↓ →       89.8333     89.6667     89.5  …  -89.6667  -89.8333  -90.0
 -180.0    -9999.0     -9999.0     -9999.0     -15.399   -13.805   -14.046
 -179.833  -9999.0     -9999.0     -9999.0     -15.9605  -14.607   -14.5545
    ⋮                                       ⋱                        ⋮
  179.5    -9999.0     -9999.0     -9999.0     -18.2955  -16.7583  -16.72
  179.667  -9999.0     -9999.0     -9999.0     -18.2847  -16.7513  -16.72
  179.833  -9999.0     -9999.0     -9999.0  …  -17.1478  -15.4243  -15.701
```


### set {#set}

`set` can be used to modify the nested properties of an objects dimensions, that are more difficult to change with `rebuild`. `set` works on the principle that dimension properties can only be in one specific field, so we generally don&#39;t have to specify which one it is. `set` will also try to update anything affected by a change you make.

This will set the `X` axis to specify points, instead of intervals:

```julia
using Rasters: Points
set(A, X => Points)
```


```
┌ 2160×1080 Raster{Union{Missing, Float32}, 2} bio5 ┐
├───────────────────────────────────────────────────┴──────────────────── dims ┐
  ↓ X Projected{Float64} -180.0:0.16666666666666666:179.83333333333331 ForwardOrdered Regular Points,
  → Y Projected{Float64} 89.83333333333333:-0.16666666666666666:-90.0 ReverseOrdered Regular Intervals{Start}
├──────────────────────────────────────────────────────────────────── metadata ┤
  Metadata{Rasters.GDALsource} of Dict{String, Any} with 1 entry:
  "filepath" => "./WorldClim/BioClim/wc2.1_10m_bio_5.tif"
├────────────────────────────────────────────────────────────────────── raster ┤
  missingval: missing
  extent: Extent(X = (-180.0, 179.83333333333331), Y = (-90.0, 90.0))
  crs: GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.25722...
└──────────────────────────────────────────────────────────────────────────────┘
    ↓ →    89.8333    89.6667    89.5       …  -89.6667  -89.8333  -90.0
 -180.0      missing    missing    missing     -15.399   -13.805   -14.046
 -179.833    missing    missing    missing     -15.9605  -14.607   -14.5545
    ⋮                                       ⋱                        ⋮
  179.5      missing    missing    missing     -18.2955  -16.7583  -16.72
  179.667    missing    missing    missing     -18.2847  -16.7513  -16.72
  179.833    missing    missing    missing  …  -17.1478  -15.4243  -15.701
```


We can also reassign dimensions, here `X` becomes `Z`:

```julia
set(A, X => Z)
```


```
┌ 2160×1080 Raster{Union{Missing, Float32}, 2} bio5 ┐
├───────────────────────────────────────────────────┴──────────────────── dims ┐
  ↓ Z Projected{Float64} -180.0:0.16666666666666666:179.83333333333331 ForwardOrdered Regular Intervals{Start},
  → Y Projected{Float64} 89.83333333333333:-0.16666666666666666:-90.0 ReverseOrdered Regular Intervals{Start}
├──────────────────────────────────────────────────────────────────── metadata ┤
  Metadata{Rasters.GDALsource} of Dict{String, Any} with 1 entry:
  "filepath" => "./WorldClim/BioClim/wc2.1_10m_bio_5.tif"
├────────────────────────────────────────────────────────────────────── raster ┤
  missingval: missing
  extent: Extent(Z = (-180.0, 179.99999999999997), Y = (-90.0, 90.0))
  crs: GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.25722...
└──────────────────────────────────────────────────────────────────────────────┘
    ↓ →    89.8333    89.6667    89.5       …  -89.6667  -89.8333  -90.0
 -180.0      missing    missing    missing     -15.399   -13.805   -14.046
 -179.833    missing    missing    missing     -15.9605  -14.607   -14.5545
    ⋮                                       ⋱                        ⋮
  179.5      missing    missing    missing     -18.2955  -16.7583  -16.72
  179.667    missing    missing    missing     -18.2847  -16.7513  -16.72
  179.833    missing    missing    missing  …  -17.1478  -15.4243  -15.701
```


`setcrs(A, crs)` and `setmappedcrs(A, crs)` will set the crs value/s of an object to any `GeoFormat` from GeoFormatTypes.jl.
