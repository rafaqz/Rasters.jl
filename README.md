# GeoData

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://rafaqz.github.io/GeoData.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://rafaqz.github.io/GeoData.jl/dev)
[![Build Status](https://travis-ci.org/rafaqz/GeoData.jl.svg?branch=master)](https://travis-ci.org/rafaqz/GeoData.jl)
[![Codecov](https://codecov.io/gh/rafaqz/GeoData.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/rafaqz/GeoData.jl)
[![Quality Assurance](https://img.shields.io/badge/GeoData.jl-%F0%9F%8C%A2-aqua.svg)](https://github.com/rafaqz/GeoData.jl)

GeoData.jl defines common types and methods for working with spatial data,
such as 2 or multidimensional raster arrays, multi-array stacks, and series of
stacks or arrays spread over multiple files. It provides a standardised
interface that allows many source data types to be used with identical syntax.

GeoData.jl extends
[DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl) so that
spatial data can be indexed using named dimensions like `Lat` and `Lon`, `Ti`
(time), which can also be used in most `Base` and `Statistics` methods like
`mean` and `reduce` where `dims` arguments are required. Much of the behaviour
is covered in the [DimensionalData
docs](https://rafaqz.github.io/DimensionalData.jl/stable/).

It provides general types `GeoArray`, `GeoStack`, and `GeoSeries`. 

R `.grd` files can be loaded natively using `GRDarray` and `GRDstack`. 
GDAL files can be loaded when
[ArchGDAL.jl](https://github.com/yeesian/ArchGDAL.jl) (v0.5 or higher) is
present, with ` GDALarray` and GDALstack. NetCDF similarly can be loaded when
[NCDatasets.jl](https://github.com/Alexander-Barth/NCDatasets.jl) is loaded,
with `NCDarray` and `NCDstack`.

When HDF5.jl is loaded, files from the Soil Moisture Active Passive
([SMAP](https://smap.jpl.nasa.gov/)) dataset can be loaded using `SMAPstack`
or `SMAPseries` to load whole directories. This is both useful for users of
SMAP, and a demonstration of the potential to build standardised interfaces 
for custom spatial dataset formats like those used in SMAP.

Files can be written to disk in all formats using `write`.

## Examples

We'll load a file from disk, and do some manipulations and plotting.

Load GeoData, and NCDatasets, download file and load it to 
an array. This netcdf file only has one layer, if it has more we 
could use `NCDstack` instead.

```julia
julia> using GeoData, NCDatasets

julia> url = "https://www.unidata.ucar.edu/software/netcdf/examples/tos_O1_2001-2002.nc";

julia> filename = download(url, "tos_O1_2001-2002.nc");

julia> A = NCDarray(filename)
NCDarray (named tos) with dimensions:
 Longitude (type Lon): Float64[1.0, 3.0, …, 357.0, 359.0] (Converted: Ordered Regular Intervals)
 Latitude (type Lat): Float64[-79.5, -78.5, …, 88.5, 89.5] (Converted: Ordered Regular Intervals)
 Time (type Ti): DateTime360Day[DateTime360Day(2001-01-16T00:00:00), DateTime360Day(2001-02-16T00:00:00), …, DateTime360Day(2002-11-16T00:00:00), DateTime360Day(2002-12-16T00:00:00)] (Sampled: Ordered Irregular Intervals)
and data: 180×170×24 Array{Union{Missing, Float32},3}
[:, :, 1]
 missing  missing     missing  …  271.437  271.445  271.459
 missing  missing     missing     271.438  271.445  271.459
...
```

Now plot every third month in the first year, just using the regular index:

```julia
julia> using Plots

juila> A[Ti(1:3:12)] |> plot
```

![Global ocean surface temperatures](https://raw.githubusercontent.com/rafaqz/GeoData.jl/media/four_pane_map.png)

Now plot Australia in the first month of 2001. Notice we are using tat/ton coordinates 
and date/time instead of regular indexes:

```julia
julia> A[Ti(Contains(DateTime360Day(2001, 01, 17))), 
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
 Time (type Ti): DateTime360Day[2001-01-16T00:00:00] (Sampled: Ordered Irregular Intervals)
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
the object doesn't have both `Lat` and `Lon` dimensions. So (as a random example) we 
could plot a transect of ocean surface temperature at 20 degree latitude :

```julia
A[Lat(Contains(20.0)), Ti(1)] |> plot
```

GeoData.jl provides a range of other methods that are being added to over time.
One example is `aggregate`, that can aggregate `GeoArray` by axis-specific amounts:

```julia
julia> aggregate(mean, A, (Ti(12), Lat(20), Lon(20))

GeoArray (named tos) with dimensions:
 Longitude (type Lon): Float64[21.0, 61.0, …, 301.0, 341.0] (Converted: Ordered Regular Intervals)
 Latitude (type Lat): Float64[-69.5, -49.5, …, 50.5, 70.5] (Converted: Ordered Regular Intervals)
 Time (type Ti): DateTime360Day[2001-01-16T00:00:00, 2002-01-16T00:00:00] (Sampled: Ordered Irregular Intervals)
and data: 9×8×2 Array{Union{Missing, Float32},3}
[:, :, 1]
 missing  277.139        missing     missing     missing     missing  missing  missing
 missing  277.126        missing     missing     missing     missing  missing  missing
```

This will also work for entire `GeoStacks` and `GeoSeries` using the same syntax.



![Temperatures at lattitude 20-21](https://raw.githubusercontent.com/rafaqz/GeoData.jl/media/lat_20.png)


## Works in progress

- Integration with Vector/DataFrame spatial types and point/line/polygon data
  types. It should be possible to select polygons of data, and convert between
  linear datasets and array formats.
- Standardised handling and conversion of spatial metadata between data formats
- Handling complex projections: Affine transformation of dimensions to indices.
  AffineMaps will be stored as a wrapper dimension in `dims`.
- Load and write the NetCDF projection format.

