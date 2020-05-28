using DimensionalData: refdims_title

@recipe function f(A::GeoArray{T,3,<:Tuple{<:Lat,<:Lon,D}}) where {T,D}
    A = prepare(A)
    nplots = size(A, 3)
    if nplots > 1
        :layout --> nplots
        :xlabel --> permutedims(string.(val(dims(A, D))))
        for i in 1:nplots
            @series begin
                seriestype := :heatmap
                aspect_ratio := 1
                subplot := i
                slice = A[:, :, i]
                lat, lon = map(maybe_reproject, dims(A))
                lon, lat, data(slice)
            end
        end
    else
        A[:, :, 1]
    end
end

# TODO generalise for any dimension order
@recipe function f(A::GeoArray{T,3}) where {T,D}
    if all(hasdim(A, (Lat, Lon)))
        permutedims(A, (Lat, Lon, Dimension))
    else
        error("Cannot plot 3 dimensional GeoArray without both Lat/Lon dims")
    end
end

@recipe function f(A::GeoArray{T,2,<:Tuple{<:Lat,<:Lon}}) where T
    A = prepare(A)
    :seriestype --> :heatmap
    :aspect_ratio --> 1
    :grid --> false
    :colorbar_title --> name(A)
    :title --> refdims_title(A)
    lat, lon = map(maybe_reproject, dims(A, (Lat, Lon)))
    lon, lat, data(A)
end

@recipe function f(A::GeoArray{T,2,<:Tuple{<:Lon,<:Lat}}) where T
    permutedims(A)
end

@recipe function f(A::AbstractGeoArray)
    GeoArray(A)
end

maybe_reproject(dim::Dimension) = maybe_reproject(mode(dim), dim)
maybe_reproject(mode::IndexMode, dim::Dimension) = val(dim)
maybe_reproject(mode::Projected, dim::Dimension) =
    maybe_reproject(crs(mode), usercrs(mode), dim)
maybe_reproject(crs, usercrs, dim::Dimension) = val(dim)
maybe_reproject(crs::GeoFormat, usercrs::GeoFormat, dim::Dimension) =
    reproject(crs, usercrs, dim, val(dim))

prepare(A) = A |> forwardorder |> maybenanmissing

maybenanmissing(A::AbstractArray{<:AbstractFloat}) = replace_missing(A, missing)
maybenanmissing(A) = A
