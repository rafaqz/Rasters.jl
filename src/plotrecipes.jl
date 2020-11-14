using DimensionalData: refdims_title

struct GeoPlot end

# We only look at arrays with Lat/Lon here.
# Otherwise they fall back to DimensionalData.jl recipes
@recipe function f(A::AbstractGeoArray)
    A = GeoArray(A)
    if all(hasdim(A, (Lat(), Lon())))
        # Heatmap or multiple heatmaps. Use GD recipes.
        GeoPlot(), prepare(A)
    else
        # This is not a Lat/Lon heatmap. Fall back to DD recipes after reprojecting
        da = A |> GeoArray |> a -> DimArray(a; dims=_maybe_mapped(dims(a)))
        DimensionalData.DimensionalPlot(), da
    end
end

# Plot 3d arrays as multiple tiled plots
@recipe function f(::GeoPlot, A::GeoArray{T,3,<:Tuple{<:Lat,<:Lon,D}}) where {T,D}
    nplots = size(A, 3)
    if nplots > 1
        :layout --> nplots
        :xguide --> permutedims(string.(val(dims(A, D))))
        for i in 1:nplots
            @series begin
                seriestype := :heatmap
                aspect_ratio := 1
                subplot := i
                slice = A[:, :, i]
                lats, lons = map(prepare, dims(slice))
                lons, lats, parent(slice)
            end
        end
    else
        GeoPlot(), A[:, :, 1]
    end
end

# # Permute for correct Lat/Lon order
@recipe function f(::GeoPlot, A::GeoArray{T,3}) where {T}
    GeoPlot(), permutedims(A, (Lat, Lon, Dimension))
end

# # Plot a sinlge 2d map
@recipe function f(::GeoPlot, A::GeoArray{T,2,<:Tuple{<:Lat,<:Lon}}) where T
    :seriestype --> :heatmap
    :aspect_ratio --> 1
    :colorbar_title --> name(A)
    :title --> refdims_title(A)
    lats, lons = map(prepare, dims(A))
    lons, lats, parent(A)
end

# # Permute for correct Lat/Lon order
@recipe function f(::GeoPlot, A::GeoArray{T,2,<:Tuple{<:Lon,<:Lat}}) where T
    GeoPlot(), permutedims(A)
end

# Plots heatmaps pixels are centered. 
# So we should center, and use the projected value.
prepare(d::Dimension) = shiftindexloci(Center(), d) |> _maybe_mapped |> index 

# Convert arrays to a consistent missing value and Forward array order
prepare(A::AbstractGeoArray) = A |> _maybe_replace_missing |> forwardorder

_maybe_replace_missing(A::AbstractArray{<:AbstractFloat}) = replace_missing(A, NaN)
_maybe_replace_missing(A) = A

_maybe_mapped(dims::Tuple) = map(_maybe_mapped, dims)
_maybe_mapped(dim::Dimension) = _maybe_mapped(mode(dim), dim)
_maybe_mapped(mode::IndexMode, dim::Dimension) = dim
_maybe_mapped(mode::Projected, dim::Dimension) = _maybe_mapped(mappedcrs(mode), dim)
_maybe_mapped(::Nothing, dim::Dimension) = dim
_maybe_mapped(::GeoFormat, dim::Dimension) = convertmode(Mapped, dim)
