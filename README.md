# GeoData

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://rafaqz.github.io/GeoData.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://rafaqz.github.io/GeoData.jl/dev)
[![Build Status](https://travis-ci.com/rafaqz/GeoData.jl.svg?branch=master)](https://travis-ci.com/rafaqz/GeoData.jl)
[![Codecov](https://codecov.io/gh/rafaqz/GeoData.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/rafaqz/GeoData.jl)
[![Aqua.jl Quality Assurance](https://img.shields.io/badge/Aquajl-%F0%9F%8C%A2-aqua.svg)](https://github.com/JuliaTesting/Aqua.jl)

GeoData.jl defines common types and methods for reading, writing and
manipulating spatial data, currently raster arrays like GeoTIFF and NetCDF,
multi-array stacks, and series of stacks or geoarrays spread over multiple files.
It provides a standardised interface that allows many source data types to be
used with identical syntax.

Data loaded with GeoData.jl has some special properties:

- Plots are always oriented the right way.
- Data is loaded lazily wherever possible using DiskArrays.jl under the hood.
  Indexing a `GeoStack` by name is always lazy, while `view` of a `GeoArray` is
  lazy and `getindex` will load to memory. `read` can be used on any 
  object to ensure that all data is loaded to memory.
- Broadcast over disk-based objects is lazy - it will only run when the array is indexed.
    Always prefer broadcasts to explicit loops - these can be very slow with
    disk-based data.
- Regions and points selected with `Between` and `Contains` select the right
  point or whole interval no matter the order of the index or it's position in
  the cell.
- For `Projected` mode you can index in any projection you want to by setting
  the `mappedcrs` keyword on construction. You don't even need to know the
  underlying projection, the conversion is handled automatically. This means
  lat/lon `EPSG(4326)` can be used across all sources seamlessly if you need
  that.
- Packages building on GeoData.jl can treat `AbstractGeoSeries`,
  `AbstractGeoStack`, and `AbstrackGeoArray` as black boxes:
  - The data could hold tiff or netcdf files, `Array`s in memory or `CuArray`s
    on the GPU - they will all behave in the same way.
  - `AbstractGeoStack` can be a Netcdf or HDF5 file, or a `NamedTuple` of
    `GDALarray` holding `.tif` files, or all `GeoArray` in memory, but be
    treated as if they are all the same thing.
  - Modelling packages do not have to deal with the specifics of spatial file
    types directly.

GeoData.jl extends
[DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl) so that
spatial data can be indexed using named dimensions like `X`, `Y` and `Ti`
(time), which can also be used in most `Base` and `Statistics` methods like
`mean` and `reduce` where `dims` arguments are required. Much of the behaviour
is covered in the [DimensionalData
docs](https://rafaqz.github.io/DimensionalData.jl/stable/).



# Objects

GeoData.jl provides general types for holding spatial data: `GeoArray`,
`GeoStack`, and `GeoSeries`. 

- `GeoArray` acts like an array, but has attached spatial metadata and can be
  indexed with e.g. it's lattitude/longitude coordinates.
- `GeoStack` is a stack of `GeoArray` layers that share (at least some)
  dimensions. Indexing with `Symbol` names will return a `GeoArray`, indexing
  with other values will return a new `GeoStack` for the selected region, or a
  `NamedTuple` of single values.
- `GeoSeries` is a `GeoArray` of `GeoArray` or `GeoStack`. Most often a
  time-series, but can be multi-dimensional as well.


# Backends

GeoData.jl relies on ArchGDAL, NCDatasets, HDF5, and MMap for loading
disk-based files. The backend is detected automatically from the file type.

- Netcdf `.nc` files can be loaded with `GeoArray(filename)` or
  `GeoStack(filename)` when
  [NCDatasets.jl](https://github.com/Alexander-Barth/NCDatasets.jl) is imported.
- Anything GDAL can load can be loaded with `GeoArray(filename)`
  [ArchGDAL.jl](https://github.com/yeesian/ArchGDAL.jl). `GeoStack` can be made
  with a `NamedTuple` holding filename `String`s or `GeoArray`s for each stack
  layer.
- Custom HDF5 `.h5` files can also be used. Currently, files from the Soil
  Moisture Active Passive ([SMAP](https://smap.jpl.nasa.gov/)) dataset can be
  loaded with `GeoStack` or `GeoSeries`. This is both useful for users of SMAP,
  and a demonstration of the potential to build standardised interfaces for
  custom spatial formats.
- `.grd/.gri` files from R can be read and written natively using memory mapping.
  These are usually the most efficient file to work with, but are not compressed
  so have large file sizes.

Files can be written to disk in all formats using `write`, and can (with some
caveats) be written to different formats than they were loaded in, providing
file-type conversion for spatial data.


## Examples

We'll load a file from disk, and do some manipulations and plotting.

Load GeoData, and NCDatasets, download file and load it to 
an array. This netcdf file only has one layer, if it has more we 
could use `NCDstack` instead.

```julia
julia> using GeoData

julia> url = "https://www.unidata.ucar.edu/software/netcdf/examples/tos_O1_2001-2002.nc";

julia> filename = download(url, "tos_O1_2001-2002.nc");

julia> A = geoarray(filename)
...
```

Now plot every third month in the first year, just using the regular index:

```julia
julia> using Plots

juila> A[Ti(1:3:12)] |> plot
```

![Global ocean surface temperatures](https://raw.githubusercontent.com/rafaqz/GeoData.jl/media/four_pane_map.png)

Now plot Australia in the first month of 2001. Notice we are using lat/lon coordinates 
and date/time instead of regular indexes:

```julia
julia> A[Ti(Near(DateTime360Day(2001, 01, 17))), 
         Lat(Between(0.0, -50.0)), 
         Lon(Between(100.0, 160.0))] |> plot
```

![Australia regional ocean surface temperature](https://raw.githubusercontent.com/rafaqz/GeoData.jl/media/aus.png)

Now get the mean over the timespan, then save it to disk, and plot it :

```julia
julia> using Statistics
julia> mean_tos = mean(A; dims=Ti)
GeoArray (named tos) with dimensions:
 Longitude (type Lon): Float64[1.0, 3.0, …, 357.0, 359.0] (Converted: Ordered Regular Intervals)
 Latitude (type Lat): Float64[-79.5, -78.5, …, 88.5, 89.5] (Converted: Ordered Regular Intervals)
 Time (type Ti): DateTime360Day[2001-01-16T00:00:00] (Sampled: Ordered Irregular Points)
and data: 180×170×1 Array{Union{Missing, Float32},3}
[:, :, 1]
 missing  missing     missing     missing  …  271.434  271.443  271.454
 missing  missing     missing     missing     271.434  271.443  271.454
...

julia> write("mean.ncd", NCDarray, mean_tos)
    Writing netcdf...
        key: "longitude" of type: Float64
        key: "latitude" of type: Float64
        key: "time" of type: DateTime360Day
        key: "tos" of type: Union{Missing, Float32}
"mean.ncd"

julia> plot(mean_tos; color=:viridis) 
```

![Mean temperatures](https://raw.githubusercontent.com/rafaqz/GeoData.jl/media/mean.png)

Plotting recipes in DimensionalData.jl are the fallback for GedData.jl when 
the object doesn't have both `X` and `Y` dimensions. So (as a random example) we 
could plot a transect of ocean surface temperature at 20 degree latitude :

```julia
A[Y(Near(20.0)), Ti(1)] |> plot
```

![Temperatures at lattitude 20-21](https://raw.githubusercontent.com/rafaqz/GeoData.jl/media/lat_20.png)


GeoData.jl provides a range of other methods that are being added to over time.
Where applicable these methods read and write lazily to and from disk-based
arrays of common raster file types. These methods also work for entire
`GeoStacks` and `GeoSeries` using the same syntax.

- `agregate` and `aggregate!` aggregate data by the same or different amounts for each axis.
- `disaggregate` and `disaggregate!` similarly disaggregate data.
- `resample` can resample data to a different size and projection, and snap to
    an existing `AbstractGeoArray`, using `warp` to access `gdalwarp`.
- `mask` and `mask!` mask and object by a polygon or GeoArray along `X/Y`, or
    arbitrary, dimensions.
- `classify` and  `classify!` classify values into categories.
- `mosaic` and `mosaic!` join rasters covering different extents into a single
    array or file.
- `crop` and `extend` shrink or extend objects to specific dimension sizes or
    the exten of another object.
- `trim` trims areas of missing values across arbitrary dimensions and stack layers.

For example, `aggregate`:

```julia
using GeoData
A = GeoArray(WorldClim{Climate}, :prec; month=1)
aggregate(mean, A, (Ti(12), Y(20), X(20))
```
