
reorderdims(dims) = map(d -> indexorder(d) == Reverse() ? rebuild(d, reverse(val(d))) : d, dims)

preparedata(A) = A |> forwardorder |> maybenanmissing |> data

maybenanmissing(A::AbstractArray{<:AbstractFloat}) = replace_missing(A, missing)
maybenanmissing(A) = A

forwardorder(A) = begin
    for (i, dim) in enumerate(dims(A))
        if arrayorder(dim) == Reverse()
            A = reverse(A; dims=dim)
        end
    end
    A
end

@recipe function f(A::AbstractGeoArray)
    GeoArray(A)
end

@recipe function f(A::GeoArray{T,3,<:Tuple{<:Lat,<:Lon,D}}) where {T,D}
    nplots = size(A, 3)
    if nplots > 1
        layout --> nplots
        xlabel --> permutedims(string.(val(dims(A, D))))
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
    seriestype --> :heatmap
    aspect_ratio --> 1
    grid --> false
    ylabel --> name(dims(A)[1])
    xlabel --> name(dims(A)[2])
    colorbar_title --> name(A)
    title --> join(map(d -> string(name(d), " ", val(d)), refdims(A)), ", ")
    (reverse(val.(reorderdims(dims(A))))..., preparedata(A))
end

@recipe function f(A::GeoArray{T,2,<:Tuple{<:Lon,<:Lat}}) where T
    permutedims(A)
end
