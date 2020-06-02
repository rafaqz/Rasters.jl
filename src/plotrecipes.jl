using DimensionalData: refdims_title

struct GeoPlot end


# We only look at arrays with Lat/Lon here.
# Otherwise they fall back to DimensionalData.jl recipes
@recipe function f(A::AbstractGeoArray)
    A = GeoArray(A)
    if all(hasdim(A, (Lat(), Lon())))
        GeoPlot(), A
    else
        println("DimPlot")
        da = A |> GeoArray |> a -> DimensionalArray(a; dims=maybe_reproject(dims(a)))
        DimensionalData.DimensionalPlot(), da
    end
end

# Plot 3d arrays as multiple tiled plots
@recipe function f(::GeoPlot, A::GeoArray{T,3,<:Tuple{<:Lat,<:Lon,D}}) where {T,D}
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
                lat, lon = maybe_reproject(dims(A))
                val(lon), val(lat), data(slice)
            end
        end
    else
        println("1")
        GeoPlot(), A[:, :, 1]
    end
end

# # Permute for correct Lat/Lon order
@recipe function f(::GeoPlot, A::GeoArray{T,3}) where {T}
    GeoPlot(), permutedims(A, (Lat, Lon, Dimension))
end

# # Plot a sinlge 2d map
@recipe function f(::GeoPlot, A::GeoArray{T,2,<:Tuple{<:Lat,<:Lon}}) where T
    A = prepare(A)
    :seriestype --> :heatmap
    :aspect_ratio --> 1
    :colorbar_title --> name(A)
    :title --> refdims_title(A)
    lat, lon = maybe_reproject(dims(A))
    val(lon), val(lat), data(A)
end

# # Permute for correct Lat/Lon order
@recipe function f(::GeoPlot, A::GeoArray{T,2,<:Tuple{<:Lon,<:Lat}}) where T
    GeoPlot(), permutedims(A)
end

maybe_reproject(dims::Tuple) = map(maybe_reproject, dims)
maybe_reproject(dim::Dimension) = maybe_reproject(mode(dim), dim)
maybe_reproject(mode::IndexMode, dim::Dimension) = dim
maybe_reproject(mode::Projected, dim::Dimension) =
    maybe_reproject(crs(mode), usercrs(mode), dim)
maybe_reproject(crs, usercrs, dim::Dimension) = dim
maybe_reproject(crs::GeoFormat, usercrs::GeoFormat, dim::Dimension) = begin
    rebuild(dim, reproject(crs, usercrs, dim, val(dim)))
end

prepare(A) = A |> forwardorder |> maybenanmissing

maybenanmissing(A::AbstractArray{<:AbstractFloat}) = replace_missing(A, missing)
maybenanmissing(A) = A
