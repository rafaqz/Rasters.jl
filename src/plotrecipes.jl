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
@recipe function f(A::GeoArray{T,3,<:Tuple{Vararg{Union{<:Lon,<:Lat,D}}}}) where {T,D}
    permutedims(A, (Lat, Lon, D))
end

@recipe function f(A::GeoArray{T,2,<:Tuple{<:Lat,<:Lon}}) where T
    A = prepare(A)
    :seriestype --> :heatmap
    :aspect_ratio --> 1
    :grid --> false
    :colorbar_title --> name(A)
    :title --> refdims_title(A)
    println(val.(dims(A, (Lat, Lon))))
    lat, lon = map(maybe_reproject, dims(A, (Lat, Lon)))
    lon, lat, data(A)
end

@recipe function f(A::AbstractGeoArray)
    GeoArray(A)
end

maybe_reproject(dim::Dimension) = maybe_reproject(mode(dim), dim, val(dim))
maybe_reproject(mode::IndexMode, dim::Dimension, vals::AbstractArray) = vals
maybe_reproject(mode::ProjectedIndex, dim::Dimension, vals::AbstractArray) =
    maybe_reproject(crs(mode), usercrs(mode), dim, vals)
maybe_reproject(crs, usercrs, dim::Dimension, vals::AbstractArray) = vals
maybe_reproject(crs::GeoFormat, usercrs::GeoFormat, dim::Dimension, vals::AbstractArray) =
    reproject(crs, usercrs, dim, vals)

@recipe function f(A::GeoArray{T,2,<:Tuple{<:Lon,<:Lat}}) where T
    permutedims(A)
end

prepare(A) = A |> forwardorder |> maybenanmissing

maybenanmissing(A::AbstractArray{<:AbstractFloat}) = replace_missing(A, missing)
maybenanmissing(A) = A
