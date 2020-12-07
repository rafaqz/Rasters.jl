struct GeoPlot end

const GeoDim = Union{GeoXDim,GeoYDim,GeoZDim}

function dim_aware_permutedims(A)
    if all(hasdim(A, (GeoXDim, GeoYDim)))
       return permutedims(A, (GeoXDim, GeoYDim))
    elseif all(hasdim(A, (GeoXDim, GeoZDim)))
        return permutedims(A, (GeoXDim, GeoZDim))
    elseif all(hasdim(A, (GeoYDim, GeoZDim)))
        return permutedims(A, (GeoYDim, GeoZDim))
    else
        return permutedims(A, (GeoXDim, GeoYDim))
    end
end

# We only look at arrays with <:GeoXDim/<:GeoYDim here.
# Otherwise they fall back to DimensionalData.jl recipes
@recipe function f(A::AbstractGeoArray)
    A = GeoArray(A)
    if all(hasdim(A, (GeoXDim, GeoYDim))) || all(hasdim(A, (GeoXDim, GeoZDim))) || all(hasdim(A, (GeoYDim, GeoZDim))) || hasdim(A, GeoZDim)
        # Heatmap or multiple heatmaps. Use GD recipes.
        GeoPlot(), prepare(A)
    else
        # This is not a GeoXDim/GeoYDim heatmap. Fall back to DD recipes after reprojecting
        da = A |> GeoArray |> a -> DimArray(a; dims=_maybe_mapped(dims(a)))
        DimensionalData.DimensionalPlot(), da
    end
end

# Plot 3d arrays as multiple tiled plots
@recipe function f(::GeoPlot, A::GeoArray{T,3,<:Tuple{<:GeoDim,<:GeoDim,D}}) where {T,D}
    nplots = size(A, 3)
    if nplots > 1
        :layout --> nplots
        :title --> permutedims(string.(val(dims(A, D))))
        for i in 1:nplots
            @series begin
                seriestype := :heatmap
                subplot := i
                slice = A[:, :, i]
                x1, x2 = map(prepare, dims(slice))
                x1, x2, parent(A)
            end
        end
    else
        GeoPlot(), A[:, :, 1]
    end
end

# # Plot a sinlge 2d map
@recipe function f(::GeoPlot, A::GeoArray{T,2,<:Tuple{<:GeoDim,<:GeoDim}}) where T
    # If colorbar is close to symmetric (< 25% difference)
    # then use a symmetric colormap and set symmetric limits
    # so zero shows up as a neutral color.
    A_min, A_max = extrema(A)
    if (A_min + A_max) / abs(A_max - A_min) < 0.25
        A_limit = max(abs(A_min), abs(A_max))
        clims = (-A_limit, A_limit)
        :seriescolor --> :balance
    else
        clims = (A_min, A_max)
    end

    clims = get(plotattributes, :clims, clims)

    if get(plotattributes, :seriestype, :none) == :contourf
        :linewidth --> 0
        :levels --> range(clims[1], clims[2], length=20)
    end

    dim1 = dims(A, 1)
    dim2 = dims(A, 2)

    xguide = name(dim1) |> string
    yguide = name(dim2) |> string
    array_name = name(A) |> string

    if haskey(dim1.metadata, :units)
        xguide *= " ($(dim1.metadata[:units]))"
    end

    if haskey(dim2.metadata, :units)
        yguide *= " ($(dim2.metadata[:units]))"
    end

    if haskey(A.metadata, :units)
        colorbar_title = array_name * " ($(A.metadata[:units]))"
    end

    :seriestype --> :heatmap
    :title --> "$array_name: $(DD._refdims_title(A))"
    :xguide --> xguide
    :yguide --> yguide
    :clims --> clims
    :colorbar_title --> array_name

    x1, x2 = map(prepare, dims(A))

    if get(plotattributes, :seriestype, :none) == :contourf
        x1, x2, clamp.(A, clims[1], clims[2])
    else
        x1, x2, parent(A)
    end
end

# # Plot a vertical 1d line
@recipe function f(::GeoPlot, A::GeoArray{T,1,<:Tuple{<:GeoZDim}}) where T
    z_dim = dims(A, 1)
    yguide = name(z_dim) |> string
    if haskey(z_dim.metadata, :units)
        yguide *= " ($(z_dim.metadata[:units]))"
    end

    xguide = name(A) |> string
    if haskey(A.metadata, :units)
        xguide *= " ($(A.metadata[:units]))"
    end

    :title --> "$(name(A)): $(DD._refdims_title(A))"
    :xguide --> xguide
    :yguide --> yguide
    :label --> ""
    z = map(prepare, dims(A))
    parent(A), z
end

# Plots heatmaps pixels are centered.
# So we should center, and use the projected value.
prepare(d::Dimension) = shiftindexloci(Center(), d) |> _maybe_mapped |> index

# Convert arrays to a consistent missing value and Forward array order
prepare(A::AbstractGeoArray) =
    _maybe_replace_missing(A) |>
    A -> reorder(A, ForwardIndex) |>
    A -> reorder(A, ForwardRelation) |>
    A -> permutedims(A, DimensionalData.commondims((GeoXDim, GeoYDim, GeoZDim, TimeDim), dims(A)))

_maybe_replace_missing(A::AbstractArray{<:AbstractFloat}) = replace_missing(A, eltype(A)(NaN))
_maybe_replace_missing(A) = A

_maybe_mapped(dims::Tuple) = map(_maybe_mapped, dims)
_maybe_mapped(dim::Dimension) = _maybe_mapped(mode(dim), dim)
_maybe_mapped(mode::IndexMode, dim::Dimension) = dim
_maybe_mapped(mode::Projected, dim::Dimension) = _maybe_mapped(mappedcrs(mode), dim)
_maybe_mapped(::Nothing, dim::Dimension) = dim
_maybe_mapped(::GeoFormat, dim::Dimension) = convertmode(Mapped, dim)
