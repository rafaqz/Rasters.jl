export resample

"""
	resample(A::AbstractGeoArray, resolution::Number;
			 crs::GeoFormat=crs(A), method::String="near")
	resample(A::AbstractGeoArray, snap::AbstractGeoArray; method::String="near")

`resample` uses `ArchGDAL.gdalwarp` to resample an `AbstractGeoArray`.

## Arguments
- `A`: The `AbstractGeoArray` to resample.
- `resolution`: A `Number` specifying the resolution for the output. 
  If the keyword argument `crs` (described below) is specified, `resolution` must be in units of the `crs`.
- `snap`: an `AbstractGeoArray` whos resolution, crs and bounds will be snapped to. 
  For best results it should roughly cover the same extent, or a subset of `A`.

## Keyword Arguments
- `crs`: A `GeoFormatTypes.GeoFormat` specifying an output crs 
  (`A` with be reprojected to `crs` in addition to being resampled). Defaults to `crs(A)`
- `method`: A `String` specifying the method to use for resampling. Defaults to `"near"` 
  (nearest neighbor resampling). See [resampling method](https://gdal.org/programs/gdalwarp.html#cmdoption-gdalwarp-r) 
  in the gdalwarp docs for a complete list of possible values.

"""
resample

function resample(A::AbstractGeoArray, resolution::Number;
				  crs::GeoFormat=crs(A),
                  method::String="near")
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

function resample(A::AbstractGeoArray, snap::AbstractGeoArray; method::String="near")
    wkt = convert(String, convert(WellKnownText, crs(snap)))
    latres, lonres = map(abs ∘ step, span(snap, (Lat(), Lon())))
    (latmin, latmax), (lonmin, lonmax) = bounds(snap, (Lat(), Lon()))
    flags = ["-t_srs", "$(wkt)",
             "-tr", "$(latres)", "$(lonres)",
             "-te", "$(lonmin)", "$(latmin)", "$(lonmax)", "$(latmax)",
             "-r", "$(method)"]
    AG.Dataset(A) do dataset
        AG.gdalwarp([dataset], flags) do warped
            GeoArray(warped)
        end
    end
end
