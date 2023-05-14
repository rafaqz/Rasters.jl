"""
	resample(x; to, size, res, method)
    resample(xs...; to=first(xs), size, res, method)

`resample` uses `ArchGDAL.gdalwarp` to resample a [`Raster`](@ref) or
[`RasterStack`](@ref) to a new `resolution` and optionally new `crs`,
or to snap to the bounds, resolution and crs of the object `to`.

# Arguments

- `x`: the object/s to resample.

# Keywords

- `to`: a `Raster`, `RasterStack`, `Tuple` of `Dimension` or `Extents.Extent`.
    If no `to` object is provided the extent will be calculated from `x`,
$RES_KEYWORD
$SIZE_KEYWORD
- `crs`: A `GeoFormatTypes.GeoFormat` coordinate reference system for the output raster, 
    such as `EPSG(x)` or `WellKnownText(string)`. Defaults to `crs(A)`.
- `method`: A `Symbol` or `String` specifying the method to use for resampling.
    From the docs for [`gdalwarp`](https://gdal.org/programs/gdalwarp.html#cmdoption-gdalwarp-r):
    * `:near`: nearest neighbour resampling (default, fastest algorithm, worst interpolation quality).
    * `:bilinear`: bilinear resampling.
    * `:cubic`: cubic resampling.
    * `:cubicspline`: cubic spline resampling.
    * `:lanczos`: Lanczos windowed sinc resampling.
    * `:average`: average resampling, computes the weighted average of all non-NODATA contributing pixels.
        rms root mean square / quadratic mean of all non-NODATA contributing pixels (GDAL >= 3.3)
    * `:mode`: mode resampling, selects the value which appears most often of all the sampled points.
    * `:max`: maximum resampling, selects the maximum value from all non-NODATA contributing pixels.
    * `:min`: minimum resampling, selects the minimum value from all non-NODATA contributing pixels.
    * `:med`: median resampling, selects the median value of all non-NODATA contributing pixels.
    * `:q1`: first quartile resampling, selects the first quartile value of all non-NODATA contributing pixels.
    * `:q3`: third quartile resampling, selects the third quartile value of all non-NODATA contributing pixels.
    * `:sum`: compute the weighted sum of all non-NODATA contributing pixels (since GDAL 3.1)

    Where NODATA values are set to `missingval`.

Note:
- GDAL may cause some unexpected changes in the data, such as returning a reversed Y dimension or
  changing the `crs` type from `EPSG` to `WellKnownText` (it will represent the same CRS).

# Example

Resample a WorldClim layer to match an EarthEnv layer:

```jldoctest
using Rasters, Plots
A = Raster(WorldClim{Climate}, :prec; month=1)
B = Raster(EarthEnv{HabitatHeterogeneity}, :evenness)

a = plot(A)
b = plot(resample(A; to=B))

savefig(a, "build/resample_example_before.png");
savefig(b, "build/resample_example_after.png"); nothing
# output
```

### Before `resample`:

![before resample](resample_example_before.png)

### After `resample`:

![after resample](resample_example_after.png)

$EXPERIMENTAL
"""
function resample end
resample(x, res; kw...) = resample(x; res, kw...)
resample(xs::RasterStackOrArray...; kw...) = resample(xs; kw...)
function resample(ser::AbstractRasterSeries, args...; kw...)
    map(x -> resample(x, args...; kw...), ser)
end
function resample(xs::Union{Tuple,NamedTuple}; to=first(xs), kw...)
    map(x -> resample(x; to, kw...), xs)
end
function resample(x::RasterStackOrArray; 
    # We need to combine the `size` and `res` keywords with 
    # the extent in extent2dims, even if we already have dims.
    to=nothing, res=nothing, crs=nothing, size=nothing, method=:near, kw...
)
    (isnothing(size) || isnothing(res)) || _size_and_res_error()

    # Flags to send to `warp`, then to GDAL
    flags = Dict{Symbol,Any}()

    # Method
    flags[:r] = method

    # Extent
    if to isa Extents.Extent || isnothing(to) || isnothing(dims(to))
        to = isnothing(to) || to isa Extents.Extent ? to : GeoInterface.extent(to)
        if !isnothing(to)
            # Get the extent of geometries
            (xmin, xmax), (ymin, ymax) = to[(:X, :Y)]
            flags[:te] = [xmin, ymin, xmax, ymax]
        end
    else
        all(hasdim(to, (XDim, YDim))) || throw(ArgumentError("`to` mush have both XDim and YDim dimensions to resize with GDAL"))
        if sampling(to, XDim) isa Points
            to = set(to, dims(to, XDim) => Intervals(Start()))
        end
        if sampling(to, YDim) isa Points
            to = set(to, dims(to, YDim) => Intervals(Start()))
        end

        # Set res from `to` if it was not already set
        if isnothing(res) && isnothing(size)
            xres, yres = map(abs ∘ step, span(to, (XDim, YDim)))
            flags[:tr] = [yres, xres]
        end
        (xmin, xmax), (ymin, ymax) = bounds(to, (XDim, YDim))
        flags[:te] = [xmin, ymin, xmax, ymax]
    end

    # CRS
    crs = if isnothing(crs) 
        if to isa Extents.Extent
            nothing
        else
            # get crs from `to` or `x` if none was passed in
            isnothing(Rasters.crs(to)) ? Rasters.crs(x) : Rasters.crs(to)
        end
    else
        crs
    end
    if !isnothing(crs)
        wkt = convert(String, convert(WellKnownText, crs))
        flags[:t_srs] = wkt
        if isnothing(Rasters.crs(x))
            @warn "You have set a crs to resample to, but the object does not have crs so GDAL will assume it is already in the target crs. Use `newraster = setcrs(raster, somecrs)` to fix this."
        end
    end

    # Resolution
    if !isnothing(res)
        xres, yres = if res isa Real
            res, res
        elseif res isa Tuple{<:Dimension{<:Real},<:Dimension{<:Real}}
            map(val, dims(res, (YDim, XDim)))
        elseif res isa Tuple{<:Real,<:Real}
            reverse(res)
        else
            throw(ArgumentError("`res` must be a `Real`, or a 2 `Tuple` of `Real` or `Dimension`s wrapping `Real`. Got $res"))
        end
        flags[:tr] = [yres, xres]
    end

    # Size
    if !isnothing(size)
        xsize, ysize = if size isa Int
            size, size
        elseif size isa Tuple{<:Dimension{Int},<:Dimension{Int}}
            map(val, dims(size, (YDim, XDim)))
        elseif size isa Tuple{Int,Int}
            reverse(size)
        else
            throw(ArgumentError("`size` must be a `Int`, or a 2 `Tuple` of `Int` or `Dimension`s wrapping `Int`. Got $size"))
        end
        flags[:ts] = [ysize, xsize]
    end

    # resample with `warp`
    resampled = warp(x, flags; kw...)

    # Return crs to the original type, from GDAL it will always be WellKnownText
    if isnothing(crs)
        return setcrs(resampled, Rasters.crs(x))
    else
        return setcrs(resampled, crs)
    end
end
