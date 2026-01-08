##################################################################################
# Plots.jl recipes

# Method specialisation singletons.
struct RasterPlot end
struct RasterZPlot end
# We only look at arrays with X, Y, Z dims here.
# Otherwise they fall back to DimensionalData.jl recipes
@recipe function f(A::AbstractRaster)
    ddplot(A) = DimArray(A; dims=_maybe_mapped(dims(A)))
    # Resample AffineProjected
    max_res = get(plotattributes, :max_res, 1000)
    A = _maybe_resample(A)
    if !(get(plotattributes, :seriestype, :none) in (:none, :heatmap, :contourf))
        DD.DimensionalPlot(), ddplot(A)
    elseif all(hasdim(A, (SpatialDim, SpatialDim)))
        # Heatmap or multiple heatmaps. Use Rasters recipes.
        A = _prepare_plots(_subsample(A, max_res))
        RasterPlot(), A
    elseif hasdim(A, ZDim) && ndims(A) == 1
        # Z dim plot, but for spatial data we want Z on the Y axis
        RasterZPlot(), _prepare_plots(A)
    else
        # Otherwise use DD recipes
        DD.DimensionalPlot(), ddplot(A)
    end
end

# Plot a sinlge 2d map
@recipe function f(::RasterPlot, A::AbstractRaster{T,2,<:Tuple{D1,D2}}) where {T,D1<:SpatialDim,D2<:SpatialDim}
    # If colorbar is close to symmetric (< 25% difference) use a symmetric
    # colormap and set symmetric limits so zero shows up as a neutral color.

    yguide, xguide = label(dims(A))

    y, x = map(_prepare_plots, dims(A))

    rdt = DD.refdims_title(A; issingle=true)
    :title --> (rdt === "" ? _maybename(A) : _maybename(A) * "\n" * rdt)
    :xguide --> xguide
    :xrotation --> -45
    :yguide --> yguide
    :grid --> true
    :gridalpha --> 0.2
    # :guidefontsize --> 10
    # :titlefontsize --> 10
    # :tickfontsize --> 6
    # :colorbar_title --> name(A)
    :linewidth --> 1
    :colorbar_titlefontsize --> 9
    :colorbar_tickfontcolor --> RGB(0.3)
    :tickfontcolor --> RGB(0.3)
    :tickcolor --> RGB(0.3)
    :tick_direction --> :out
    :framestyle --> :box
    :foreground_color_axis --> RGB(0.3)
    :foreground_color_border --> RGB(0.3)
    :seriescolor --> :batlow

    if all(d -> lookup(d) isa Mapped, (x, y))
        :xlims --> mappedbounds(x)
        :ylims --> mappedbounds(y)
        :aspect_ratio --> :equal
    else
        :xlims --> bounds(A, x)
        :ylims --> bounds(A, y)
        bnds = bounds(A, (D1, D2))
        # TODO: Rethink this....
        s1, s2 = map(((l, u),) -> (u - l), bnds) ./ (size(A, D1), size(A, D2))
        square_pixels = s2 / s1
        :aspect_ratio --> square_pixels
    end


    if eltype(A) <: ColorTypes.Colorant
        index(x), index(y), parent(A)
    elseif get(plotattributes, :seriestype, :none) == :contourf
        A = replace_missing(A, missing)
        clims = extrema(skipmissing(A))
        :levels --> range(clims[1], clims[2], length=20)
        index(x), index(y), clamp.(A, clims[1], clims[2])
    else
        :seriestype --> :heatmap
        A = replace_missing(A, missing)
        index(x), index(y), parent(A)
    end
end

# Plot a vertical 1d line
@recipe function f(::RasterZPlot, A::Raster{T,1}) where T
    z_dim = dims(A, ZDim)
    yguide = label(z_dim)
    xguide = label(A)
    rdt = DD.refdims_title(A; issingle=true)
    :title --> (rdt == "" ? _maybename(A) : _maybename(A) * "\n" * rdt)
    :xguide --> xguide
    :yguide --> yguide
    :label --> ""
    z = map(_prepare_plots, dims(A))
    parent(A), index(z)
end

# Plot 3d arrays as multiple tiled plots
@recipe function f(::RasterPlot, A::Raster{T,3,<:Tuple{<:SpatialDim,<:SpatialDim,D}}) where {T,D}
    nplots = size(A, 3)
    if nplots > 1
        :plot_title --> name(A)
        # Plot as a RasterSeries
        slice(A, dims(A, 3))
    else
        RasterPlot(), view(A, :, :, 1)
    end
end

# We only look at arrays with X, Y, Z dims here.
# Otherwise they fall back to DimensionalData.jl recipes
@recipe function f(st::AbstractRasterStack)
    D = map(d -> rebuild(d, 1), otherdims(st, (X(), Y())))
    if length(D) > 0
        st = view(st, D...)
    end
    nplots = length(keys(st))
    if nplots > 1

        if haskey(plotattributes, :layout)
            first = true
            i = 1
            for raster in values(st)
                @series begin
                    :titlefontsize := 9
                    :tickfontsize := 7
                    subplot := i
                    if !first
                        :xformatter := _ -> ""
                        :yformatter := _ -> ""
                        :xguide := ""
                        :yguide := ""
                    end
                    raster
                end
                i += 1
                first = false
            end
        else
            ncols, nrows = _balance_grid(nplots)
            :layout --> (ncols, nrows)
            rasters = values(st)
            for r in 1:nrows, c in 1:ncols
                i = (r + (c - 1) * nrows)
                @series begin
                    :titlefontsize := 9
                    :tickfontsize := 7
                    subplot := i
                    if c != ncols || r != 1
                        :xformatter := _ -> ""
                        :yformatter := _ -> ""
                        :xguide := ""
                        :yguide := ""
                    end
                    if i <= nplots
                        A = rasters[i]
                        title := string(keys(st)[i])
                        if length(dims(A, (XDim, YDim))) > 0
                            # Get a view of the first slice of the X/Y dimension
                            rasters[i]
                        else
                            framestyle := :none
                            legend := :none
                            []
                        end
                    else
                        framestyle := :none
                        legend := :none
                        []
                    end
                end
            end
        end
    else
        first(st)
    end
end

@recipe function f(A::RasterSeries{<:Any,1})
    # link --> :both
    # :colorbar := false
    if haskey(plotattributes, :layout)
        first = true
        for (i, raster) in enumerate(A)
            @series begin
                :titlefontsize := 9
                :tickfontsize := 7
                subplot := i
                if !first
                    :xformatter := _ -> ""
                    :yformatter := _ -> ""
                    :xguide := ""
                    :yguide := ""
                end
                raster
            end
            first = false
        end
    else
        thinned, plotinds, nplots = _maybe_thin_plots(A)
        titles = string.(index(A, dims(A, 1)))
        ncols, nrows = _balance_grid(nplots)
        :layout --> (ncols, nrows)
        for r in 1:nrows, c in 1:ncols
            i = (r + (c - 1) * nrows)
            @series begin
                :titlefontsize := 9
                :tickfontsize := 7
                subplot := i
                if c != ncols || r != 1
                    :xformatter := _ -> ""
                    :yformatter := _ -> ""
                    :xguide := ""
                    :yguide := ""
                end
                if i <= nplots
                    title := titles[i]
                    thinned[i]
                else
                    :framestyle := :none
                    :legend := :none
                    []
                end
            end
        end
    end
end


##################################################################################
# Makie.jl recipes
# initial definitions of `rplot`, to get around the extension package availability question

"""
    Rasters.rplot([position::GridPosition], raster; kw...)

`raster` may be a `Raster` (of 2 or 3 dimensions) or a `RasterStack` whose underlying rasters are 2 dimensional, or 3-dimensional with a singleton (length-1) third dimension.

## Keywords

- `plottype = Makie.Heatmap`: The type of plot. Can be any Makie plot type which accepts a `Raster`; in practice, `Heatmap`, `Contour`, `Contourf` and `Surface` are the best bets.
- `axistype = Makie.Axis`: The type of axis. This can be an `Axis`, `Axis3`, `LScene`, or even a `GeoAxis` from GeoMakie.jl.
- `X = XDim`: The X dimension of the raster.
- `Y = YDim`: The Y dimension of the raster.
- `Z = YDim`: The Y dimension of the raster.
- `draw_colorbar = true`: Whether to draw a colorbar for the axis or not.
- `colorbar_position = Makie.Right()`: Indicates which side of the axis the colorbar should be placed on.  Can be `Makie.Top()`, `Makie.Bottom()`, `Makie.Left()`, or `Makie.Right()`.
- `colorbar_padding = Makie.automatic`: The amount of padding between the colorbar and its axis.  If `automatic`, then this is set to the width of the colorbar.
- `title = Makie.automatic`: The titles of each plot. If `automatic`, these are set to the name of the band.
- `xlabel = Makie.automatic`: The x-label for the axis. If `automatic`, set to the dimension name of the X-dimension of the raster.
- `ylabel = Makie.automatic`: The y-label for the axis. If `automatic`, set to the dimension name of the Y-dimension of the raster.
- `colorbarlabel = ""`: Usually nothing, but here if you need it. Sets the label on the colorbar.
- `colormap = nothing`: The colormap for the heatmap. This can be set to a vector of colormaps (symbols, strings, `cgrad`s) if plotting a 3D raster or RasterStack.
- `colorrange = Makie.automatic`: The colormap for the heatmap.  This can be set to a vector of `(low, high)` if plotting a 3D raster or RasterStack.
- `nan_color = :transparent`: The color which `NaN` values should take. Default to transparent.
"""
function rplot(args...)
    @error("Please load `Makie.jl` and then call this function. If Makie is loaded, then you can't call `rplot` with no arguments!")
end

# define the theme

# this function is defined so that we can override style_rasters in RastersMakieExt
function style_rasters end
function color_rasters end

function theme_rasters()
    return merge(style_rasters(), color_rasters())
end


##################################################################################
# Utils

# Plots.jl heatmaps pixels are centered.
# So we should center the index, and use the projected value.
_prepare_plots(d::Dimension) = d |> _maybe_shift |> _maybe_mapped
# Convert arrays to a consistent missing value and Forward array order
_prepare_plots(A::AbstractRaster) = A |> _reorder |> _permute
_reorder(A) = reorder(A, DD.ForwardOrdered)
_permute(A) = permutedims(A, DD.commondims(>:, (ZDim, YDim, XDim, TimeDim, Dimension), dims(A)))

function _subsample(A, max_res)
    ssdims = dims(A, (XDim, YDim, ZDim))[1:2]
    # Aggregate based on the number of pixels
    s1, s2 = size(A, ssdims[1]), size(A, ssdims[2])
    ag = floor(Int, max(s1, s2) / max_res) + 1
    if ag == 1
        return read(A)
    else
        d1 = rebuild(ssdims[1], 1:ag:s1)
        d2 = rebuild(ssdims[2], 1:ag:s2)
        # TODO make this actually load lazily.
        # DiskArrays.jl does not handle StepRange views
        return view(read(A), d1, d2)
    end
end

function _maybe_resample(A)
    if all(hasdim(A, (X, Y)))
        _maybe_resample(first(lookup(A, (X, Y))), A)
    else
        return A
    end
end
_maybe_resample(lookup, A) = A

_maybename(A)::String = _maybename(name(A))
_maybename(::Name{N}) where N = _maybename(N)::String
_maybename(n::NoName)::String = ""
_maybename(n::Symbol)::String = String(n)
_maybename(n::AbstractString)::String = String(n)

_maybe_replace_missing(A::AbstractArray{<:AbstractFloat}) = replace_missing(A, eltype(A)(NaN))
_maybe_replace_missing(A) = A

_maybe_shift(d) = _maybe_shift(sampling(d), d)
_maybe_shift(::Intervals, d) = DD.maybeshiftlocus(Center(), d)
_maybe_shift(sampling, d) = d

_maybe_mapped(dims::Tuple) = map(_maybe_mapped, dims)
_maybe_mapped(dim::Dimension) = _maybe_mapped(lookup(dim), dim)
_maybe_mapped(lookup::Lookup, dim::Dimension) = dim
_maybe_mapped(lookup::Projected, dim::Dimension) = _maybe_mapped(mappedcrs(lookup), dim)
_maybe_mapped(::Nothing, dim::Dimension) = dim
_maybe_mapped(::GeoFormat, dim::Dimension) = convertlookup(Mapped, dim)

function _balance_grid(nplots)
    ncols = (nplots - 1) ÷ ceil(Int, sqrt(nplots)) + 1
    nrows = (nplots - 1) ÷ ncols + 1
    return ncols, nrows
end

function _maybe_thin_plots(A::AbstractRasterSeries)
    nplots = length(A)
    if nplots > 16
        plotinds = round.(Int, 1:nplots//16:nplots)
        @info "too many raster heatmaps: plotting 16 slices from $nplots"
        thinned = @views A[plotinds]
        nplots = length(plotinds)
        return thinned, plotinds, nplots
    else
        return A, collect(eachindex(A)), nplots
    end
end
