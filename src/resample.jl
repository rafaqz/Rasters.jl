export resample

"""
	resample(A::AbstractGeoArray, resolution::Number; crs, method)
	resample(A::AbstractGeoArray; to::AbstractGeoArray, method)

`resample` uses `ArchGDAL.gdalwarp` to resample an `AbstractGeoArray`.

# Arguments

- `A`: The `AbstractGeoArray` to resample.
- `resolution`: A `Number` specifying the resolution for the output.
    If the keyword argument `crs` (described below) is specified, `resolution` must be in units of the `crs`.

# Keywords

- `to`: an `AbstractGeoArray` whos resolution, crs and bounds will be snapped to.
    For best results it should roughly cover the same extent, or a subset of `A`.
- `crs`: A `GeoFormatTypes.GeoFormat` specifying an output crs
    (`A` will be reprojected to `crs` in addition to being resampled). Defaults to `crs(A)`
- `method`: A `Symbol` or `String` specifying the method to use for resampling. Defaults to `:near`
    (nearest neighbor resampling). See [resampling method](https://gdal.org/programs/gdalwarp.html#cmdoption-gdalwarp-r)
    in the gdalwarp docs for a complete list of possible values.

"""
function resample end
function resample(A::AbstractGeoArray, resolution::Number;
    crs::GeoFormat=crs(A), method=:near
)
    wkt = convert(String, convert(WellKnownText, crs))
    flags = ["-t_srs", "$(wkt)",
             "-tr", "$(resolution)", "$(resolution)",
             "-r", "$(method)"]
    AG.Dataset(A) do dataset
        AG.gdalwarp([dataset], flags) do warped
            GeoArray(warped)
        end
    end
end
function resample(A::AbstractGeoArray; to, method=:near)
    wkt = convert(String, convert(WellKnownText, crs(to)))
    latres, lonres = map(abs ∘ step, span(to, (Y(), X())))
    (latmin, latmax), (lonmin, lonmax) = bounds(to, (Y(), X()))
    flags = ["-t_srs", "$(wkt)",
             "-tr", "$(latres)", "$(lonres)",
             "-te", "$(lonmin)", "$(latmin)", "$(lonmax)", "$(latmax)",
             "-r", "$(string(method))"]
    AG.Dataset(A) do dataset
        AG.gdalwarp([dataset], flags) do warped
            read(GeoArray(warped))
        end
    end
end
