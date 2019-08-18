@recipe function f(ga::AbstractGeoArray{T,3,<:Tuple{<:Lat,<:Lon,D}}) where {T,D}
    nplots = size(ga, 3)
    if nplots > 1
        layout --> nplots
        # How to make this work?
        plot_title --> join(label(ga), label(dims(ga)[3]), " ")
        for i in 1:nplots
            @series begin
                seriestype := :heatmap
                colorbar := false
                ticks := false
                subplot := i
                replace(parent(ga[:, :, i]), missingval(ga) => NaN)
            end   
        end
    else
        ga[:, :, 1]
    end
end

@recipe function f(ga::AbstractGeoArray{T,3,<:Tuple{<:Lon,<:Lat,D}}) where {T,D}
    permutedims(ga, (Lat(), Lon(), basetype(D)()))
end

@recipe function f(ga::AbstractGeoArray{T,2,<:Tuple{<:Lat,<:Lon}}) where T
    seriestype --> :heatmap
    aspect_ratio --> 1
    grid --> false
    ylabel --> label(dims(ga)[1])
    xlabel --> label(dims(ga)[2])
    colorbar_title --> label(ga)
    title --> label(refdims(ga))
    data = replace(parent(ga), missingval(ga) => NaN)
    reverse(val.(dims(ga)))..., data
end

@recipe function f(ga::AbstractGeoArray{T,2,<:Tuple{<:Lon,<:Lat}}) where T
    permutedims(ga)
end
