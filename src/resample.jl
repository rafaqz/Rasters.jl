export resample

"""
	resample(A::AbstractGeoArray, resolution::Number;
			   crs::GeoFormat=crs(A),
			   method::String="near")

`resample` uses `ArchGDAL.gdalwarp` to resample an `AbstractGeoArray`.

## Arguments
- `A`: The `AbstractGeoArray` to resample.
- `resolution`: A `Number` specifying the resolution for the output. If the keyword argument `crs` (described below) is specified, `resolution` must be in units of the `crs`.

## Keyword Arguments
- `crs`: A `GeoFormatTypes.GeoFormat` specifying an output crs (`A` with be reprojected to `crs` in addition to being resampled). Defaults to `crs(A)`
- `method`: A `String` specifying the method to use for resampling. Defaults to `"near"` (nearest neighbor resampling). See [resampling method](https://gdal.org/programs/gdalwarp.html#cmdoption-gdalwarp-r) in the gdalwarp docs for a complete list of possible values.

"""
function resample(A::AbstractGeoArray, resolution::Number;
				  crs::GeoFormat=crs(A),
                  method::String="near")
    wkt = convert(String, convert(WellKnownText, proj))

    AG.Dataset(A) do dataset
        AG.gdalwarp([dataset], ["-t_srs", "$(wkt)",
                                "-tr", "$(resolution)", "$(resolution)",
                                "-r", "$(method)"]) do warped
            GeoArray(warped)
        end
    end
end
