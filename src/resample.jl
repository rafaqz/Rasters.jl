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
