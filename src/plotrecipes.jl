
reorderdims(dims) =
    map(dims) do d
        if isrev(DimensionalData.indexorder(d))
            rebuild(d, reverse(val(d)))
        else
            d
        end
    end

preparedata(A) = A |> maybe_nanmissing |> data

maybe_nanmissing(A::AbstractArray{<:AbstractFloat}) = replace_missing(A, missing)
maybe_nanmissing(A) = A

@recipe function f(A::AbstractGeoArray)
    GeoArray(A)
end

@recipe function f(A::GeoArray{T,3,<:Tuple{<:Lat,<:Lon,D}}) where {T,D}
    nplots = size(A, 3)
    if nplots > 1
        :layout --> nplots
        :xlabel --> permutedims(string.(val(dims(A, D))))
        for i in 1:nplots
            @series begin
                seriestype := :heatmap
                aspect_ratio := 1
                subplot := i
                (reverse(val.(reorderdims(dims(A, (Lat, Lon)))))..., preparedata(A[:, :, i]))
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
    :seriestype --> :heatmap
    :aspect_ratio --> 1
    :grid --> false
    # :colorbar_title --> name(A)
    (reverse(val.(reorderdims(dims(A))))..., preparedata(A))
    preparedata(A), DimensionalData.HeatMapLike
end
