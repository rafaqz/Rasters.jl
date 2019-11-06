
reorderdims(dims) = map(d -> indexorder(d) == Reverse() ? rebuild(d, reverse(val(d))) : d, dims)

preparedata(A) = begin
    data = parent(replace_missing(A, NaN))
    for (i, dim) in enumerate(dims(A))
        if arrayorder(dim) == Reverse()
            data = reverse(data; dims=i)
        end
    end
    data
end

@recipe function f(A::AbstractGeoArray{T,3,<:Tuple{<:Lat,<:Lon,D}}) where {T,D}
    nplots = size(A, 3)
    if nplots > 1
        layout --> nplots
        # How to make this work?
        plot_title --> join(label(A), label(dims(A)[3]), " ")
        for i in 1:nplots
            @series begin
                seriestype := :heatmap
                aspect_ratio := 1
                colorbar := false
                ticks := false
                subplot := i
                preparedata(A)
            end   
        end
    else
        A[:, :, 1]
    end
end

# TODO generalise for any dimension order
@recipe function f(A::AbstractGeoArray{T,3,<:Tuple{<:Lon,<:Lat,D}}) where {T,D}
    permutedims(A, (Lat(), Lon(), basetype(D)()))
end

@recipe function f(A::AbstractGeoArray{T,2,<:Tuple{<:Lat,<:Lon}}) where T
    seriestype --> :heatmap
    aspect_ratio --> 1
    grid --> false
    ylabel --> label(dims(A)[1])
    xlabel --> label(dims(A)[2])
    colorbar_title --> label(A)
    title --> label(refdims(A))
    reverse(val.(reorderdims(dims(A))))..., preparedata(A)
end

@recipe function f(A::AbstractGeoArray{T,2,<:Tuple{<:Lon,<:Lat}}) where T
    permutedims(A)
end
