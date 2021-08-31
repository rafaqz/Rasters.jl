# Method specialisation singletons. 
struct GeoPlot end
struct GeoZPlot end
# We only look at arrays with X, Y, Z dims here.
# Otherwise they fall back to DimensionalData.jl recipes
@recipe function f(A::AbstractGeoArray)
    ddplot(A) = DimArray(A; dims=_maybe_mapped(dims(A)))
    max_res = get(plotattributes, :max_res, 1000)
    if !(get(plotattributes, :seriestype, :none) in (:none, :heatmap, :contourf))
        DD.DimensionalPlot(), ddplot(A)
    elseif all(hasdim(A, (SpatialDim, SpatialDim)))
        # Heatmap or multiple heatmaps. Use GD recipes.
        GeoPlot(), _prepare(_subsample(A, max_res))
    elseif hasdim(A, ZDim) && ndims(A) == 1
        # Z dim plot, but for spatial data we want Z on the Y axis
        GeoZPlot(), _prepare(A)
    else
        DD.DimensionalPlot(), ddplot(A)
    end
end

# Plot 3d arrays as multiple tiled plots
@recipe function f(::GeoPlot, A::GeoArray{T,3,<:Tuple{<:SpatialDim,<:SpatialDim,D}}) where {T,D}
    nplots = size(A, 3)
    if nplots > 1
        ncols = (nplots - 1) ÷ ceil(Int, sqrt(nplots)) + 1
        nrows = (nplots - 1) ÷ ncols + 1
        :layout --> (ncols, nrows)
        # link --> :both
        # clims = extrema(A)
        colorbar := false
        titles = string.(index(A, D))
        for r in 1:nrows, c in 1:ncols
            i = (r + (c - 1) * nrows)
            @series begin
                titlefontsize := 7
                tickfontsize := 6 
                tickfontsize := 6 
                subplot := i
                if c != ncols || r != 1
                    xformatter := _ -> ""
                    yformatter := _ -> ""
                    xguide := ""
                    yguide := ""
                end
                if i <= nplots
                    title := titles[i]
                    GeoPlot(), A[:, :, i]
                else
                    framestyle := :none
                    legend := :none
                    []
                end
            end
        end
    else
        GeoPlot(), A[:, :, 1]
    end
end
# Plot a sinlge 2d map
@recipe function f(::GeoPlot, A::GeoArray{T,2,<:Tuple{<:SpatialDim,<:SpatialDim}}) where T
    # If colorbar is close to symmetric (< 25% difference) use a symmetric 
    # colormap and set symmetric limits so zero shows up as a neutral color.
    A_min, A_max = extrema(skipmissing(A))
    if (A_min + A_max) / abs(A_max - A_min) < 0.25
        A_limit = max(abs(A_min), abs(A_max))
        clims = (-A_limit, A_limit)
        :seriescolor --> :curl
    else
        clims = A_min, A_max
    end

    yguide, xguide = label(dims(A))

    rdt = DD.refdims_title(A; issingle=true)
    :title --> rdt === "" ? _maybename(A) : _maybename(A) * " " * rdt 
    :xguide --> xguide
    :yguide --> yguide
    :clims --> clims
    :axes --> :none
    :guidefontsize --> 8
    :tickfontsize --> 6 
    :titlefontsize --> 10
    :colorbar_title --> name(A)
    :colorbar_titlefontsize --> 9
    :colorbar_tickfontcolor --> RGB(0.3)
    :tickfontcolor --> RGB(0.3)
    :framestyle --> :grid
    :widen --> true
    :foreground_color_axis --> RGB(0.5)
    :seriescolor --> :magma
    :gridalpha --> 0.2

    if mappedcrs(A) === nothing
        :aspect_ratio --> :equal
    else
        bnds = bounds(A, (X, Y))
        s1, s2 = map(((l, u),) -> (u - l), bnds) ./ size(A)
        ratio = s1 / s2
        :aspect_ratio --> ratio
    end

    ys, xs = map(_prepare, dims(A))

    if get(plotattributes, :seriestype, :none) == :contourf
        :linewidth --> 0
        :levels --> range(clims[1], clims[2], length=20)
        xs, ys, clamp.(A, clims[1], clims[2])
    else
        :seriestype --> :heatmap
        xs, ys, parent(A)
    end
end

# Plot a vertical 1d line
@recipe function f(::GeoZPlot, A::GeoArray{T,1}) where T
    z_dim = dims(A, ZDim)
    yguide = label(z_dim)
    xguide = label(A)
    rdt = DD.refdims_title(A; issingle=true)
    :title --> rdt == "" ? _maybename(A) : _maybename(A) * " " * rdt 
    :xguide --> xguide
    :yguide --> yguide
    :label --> ""
    z = map(_prepare, dims(A))
    parent(A), z
end

# Plots.jl heatmaps pixels are centered.
# So we should center the index, and use the projected value.
_prepare(d::Dimension) = d |> _maybe_shift |> _maybe_mapped |> index
# Convert arrays to a consistent missing value and Forward array order
function _prepare(A::AbstractGeoArray)
    reorder(A, ForwardIndex) |> 
    a -> reorder(a, ForwardRelation) |>
    a -> permutedims(a, DD.commondims(>:, (ZDim, YDim, XDim, TimeDim, Dimension), dims(A))) |>
    a -> replace_missing(a, missing)
end

function _subsample(A, max_res)
    ssdims = dims(A, (XDim, YDim, ZDim))[1:2]
    # Aggregate based on the number of pixels
    s1, s2 = size(A, ssdims[1]), size(A, ssdims[2])
    ag = floor(Int, max(s1, s2) / max_res) + 1
    A[DD.basetypeof(ssdims[1])(1:ag:s1), DD.basetypeof(ssdims[2])(1:ag:s2)]
end

_maybename(A) = _maybename(name(A))
_maybename(n::Name{N}) where N = _maybename(N) 
_maybename(n::NoName) = ""
_maybename(n::Symbol) = string(n)

_maybe_replace_missing(A::AbstractArray{<:AbstractFloat}) = replace_missing(A, eltype(A)(NaN))
_maybe_replace_missing(A) = A

_maybe_shift(d) = _maybe_shift(sampling(d), d)
_maybe_shift(::Intervals, d) = DD.maybeshiftlocus(Center(), d)
_maybe_shift(sampling, d) = d

_maybe_mapped(dims::Tuple) = map(_maybe_mapped, dims)
_maybe_mapped(dim::Dimension) = _maybe_mapped(mode(dim), dim)
_maybe_mapped(mode::IndexMode, dim::Dimension) = dim
_maybe_mapped(mode::Projected, dim::Dimension) = _maybe_mapped(mappedcrs(mode), dim)
_maybe_mapped(::Nothing, dim::Dimension) = dim
_maybe_mapped(::GeoFormat, dim::Dimension) = convertmode(Mapped, dim)

# We don't show the Band label for a single-band raster,
# it's just not interesting information.
function DD.refdims_title(refdim::Band; issingle=false)
    if issingle
        ""
    else
        string(name(refdim), ": ", DD.refdims_title(mode(refdim), refdim))
    end
end
