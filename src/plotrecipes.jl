using DimensionalData: refdims_title

@recipe function f(A::AbstractGeoArray{T,3,<:Tuple{<:Lat,<:Lon,D}}) where {T,D}
    A = GeoArray(A)
    nplots = size(A, 3)
    if nplots > 1
        :layout --> nplots
        :xlabel --> permutedims(string.(val(dims(A, D))))
        for i in 1:nplots
            @series begin
                seriestype := :heatmap
                aspect_ratio := 1
                subplot := i
                slice = prepare(A[:, :, i])
                lat, lon = map(maybe_reproject, dims(A))
                lon, lat, data(slice)
            end
        end
    else
        A[:, :, 1]
    end
end

# TODO generalise for any dimension order
@recipe function f(A::AbstractGeoArray{T,3,<:Tuple{Vararg{Union{<:Lon,<:Lat,D}}}}) where {T,D}
    permutedims(A, (Lat, Lon, D))
end

@recipe function f(A::AbstractGeoArray{T,2,<:Tuple{<:Lat,<:Lon}}) where T
    :seriestype --> :heatmap
    :aspect_ratio --> 1
    :grid --> false
    :colorbar_title --> name(A)
    :title --> refdims_title(A)
    A = prepare(A)
    lat, lon = map(maybe_reproject, dims(A))
    lon, lat, data(A)
end

maybe_reproject(dim::Dimension) = maybe_reproject(mode(dim), dim, val(dim))
maybe_reproject(mode::IndexMode, dim::Dimension, vals::AbstractArray) = vals
maybe_reproject(mode::ProjectedIndex, dim::Dimension, vals::AbstractArray) =
    maybe_reproject(crs(mode), usercrs(mode), dim, vals)
maybe_reproject(crs, usercrs, dim::Dimension, vals::AbstractArray) = vals
maybe_reproject(crs::GeoFormat, usercrs::GeoFormat, dim::Dimension, vals::AbstractArray) =
    reproject(crs, usercrs, dim, vals)
# maybe_reproject(mode::ProjectedIndex, dim::Lat, vals::AbstractArray) =
    # [r[1] for r in ArchGDAL.reproject([(0.0, v) for v in vals], crs(mode), usercrs(mode))]
# maybe_reproject(mode::ProjectedIndex, dim::Lon, vals::AbstractArray) =
    # [r[2] for r in ArchGDAL.reproject([(v, 0.0) for v in vals], crs(mode), usercrs(mode))]

@recipe function f(A::AbstractGeoArray{T,2,<:Tuple{<:Lon,<:Lat}}) where T
    permutedims(A)
end

prepare(A) = A |> forwardorder |> maybenanmissing

maybenanmissing(A::AbstractArray{<:AbstractFloat}) = replace_missing(A, missing)
maybenanmissing(A) = A
