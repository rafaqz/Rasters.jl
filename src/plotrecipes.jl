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
                lat, lon = forwardorder(dims(A, (Lat, Lon)))
                val(lon), val(lat), preparedata(A[:, :, i])
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
    lat, lon = forwardorder(dims(A))
    val(lon), val(lat), preparedata(A)
end

@recipe function f(A::AbstractGeoArray{T,2,<:Tuple{<:Lon,<:Lat}}) where T
    permutedims(A)
end

preparedata(A) = A |> forwardorder |> maybenanmissing |> data

maybenanmissing(A::AbstractArray{<:AbstractFloat}) = replace_missing(A, missing)
maybenanmissing(A) = A

forwardorder(A::AbstractArray) = begin
    for (i, dim) in enumerate(dims(A))
        if arrayorder(dim) == Reverse()
            A = reverse(A; dims=dim)
        end
    end
    A
end

forwardorder(dims::Tuple) =
    map(dims) do d
        if isrev(DimensionalData.indexorder(d))
            rebuild(d, reverse(val(d)))
        else
            d
        end
    end
