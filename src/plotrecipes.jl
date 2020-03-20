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
                _reproject(lon), _reproject(lat), preparedata(A[:, :, i])
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
    _reproject(lon), _reproject(lat), preparedata(A)
end

_reproject(dim::Dimension) = _reproject(grid(dim), dim, val(dim))
_reproject(grid, dim::Lat, vals::AbstractArray) = 
    [r[1] for r in ArchGDAL.reproject([(0.0, v) for v in vals], crs(grid), selectorcrs(grid))]
_reproject(grid, dim::Lon, vals::AbstractArray) =                                           
    [r[2] for r in ArchGDAL.reproject([(v, 0.0) for v in vals], crs(grid), selectorcrs(grid))]

@recipe function f(A::AbstractGeoArray{T,2,<:Tuple{<:Lon,<:Lat}}) where T
    permutedims(A)
end

preparedata(A) = A |> forwardorder |> maybenanmissing |> data

maybenanmissing(A::AbstractArray{<:AbstractFloat}) = replace_missing(A, missing)
maybenanmissing(A) = A
