# Method specialisation singletons.
struct RasterPlot end
struct RasterZPlot end
# We only look at arrays with X, Y, Z dims here.
# Otherwise they fall back to DimensionalData.jl recipes
@recipe function f(A::AbstractRaster)
    ddplot(A) = DimArray(A; dims=_maybe_mapped(dims(A)))
    # Resample AffineProjected
    max_res = get(plotattributes, :max_res, 1000)
    if all(hasdim(A, (X, Y))) && first(lookup(A, (X, Y))) isa AffineProjected
        l = first(lookup(A, (X, Y)))
        res = Y(abs(l.affinemap.linear[1, 1])), X(abs(l.affinemap.linear[2, 2]))
        A = resample(A, res)
    end
    if !(get(plotattributes, :seriestype, :none) in (:none, :heatmap, :contourf))
        DD.DimensionalPlot(), ddplot(A)
    elseif all(hasdim(A, (SpatialDim, SpatialDim)))
        # Heatmap or multiple heatmaps. Use GD recipes.
        A = _prepare(_subsample(A, max_res))
        RasterPlot(), A
    elseif hasdim(A, ZDim) && ndims(A) == 1
        # Z dim plot, but for spatial data we want Z on the Y axis
        RasterZPlot(), _prepare(A)
    else
        DD.DimensionalPlot(), ddplot(A)
    end
end
# Plot a sinlge 2d map
@recipe function f(::RasterPlot, A::AbstractRaster{T,2,<:Tuple{D1,D2}}) where {T,D1<:SpatialDim,D2<:SpatialDim}
    # If colorbar is close to symmetric (< 25% difference) use a symmetric
    # colormap and set symmetric limits so zero shows up as a neutral color.

    yguide, xguide = label(dims(A))

    y, x = map(_prepare, dims(A))

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
    z = map(_prepare, dims(A))
    parent(A), index(z)
end

# Plot 3d arrays as multiple tiled plots
@recipe function f(::RasterPlot, A::Raster{T,3,<:Tuple{<:SpatialDim,<:SpatialDim,D}}) where {T,D}
    nplots = size(A, 3)
    if nplots > 1
        :plot_title --> name(A)
        ncols, nrows = _balance_grid(nplots)
        :layout --> (ncols, nrows)
        # link --> :both
        # :colorbar := false
        titles = string.(index(A, D))
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
                    RasterPlot(), A[:, :, i]
                else
                    :framestyle := :none
                    :legend := :none
                    []
                end
            end
        end
    else
        RasterPlot(), A[:, :, 1]
    end
end

function _balance_grid(nplots)
    ncols = (nplots - 1) ÷ ceil(Int, sqrt(nplots)) + 1
    nrows = (nplots - 1) ÷ ncols + 1
    return ncols, nrows
end

# We only look at arrays with X, Y, Z dims here.
# Otherwise they fall back to DimensionalData.jl recipes
@recipe function f(st::AbstractRasterStack)
    nplots = length(keys(st))
    if nplots > 1
        ncols, nrows = _balance_grid(nplots)
        :layout --> (ncols, nrows)
        # colorbar := false
        max_res = get(plotattributes, :max_res, 1000/(max(nrows, ncols)))

        l = DD.layers(st)
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
                    A = l[i]
                    title := string(keys(st)[i])
                    if length(dims(A, (XDim, YDim))) > 0
                        # Get a view of the first slice of the X/Y dimension
                        ods = otherdims(A, (X, Y))
                        if length(ods) > 0
                            od1s = map(d -> DD.basetypeof(d)(firstindex(d)), ods)
                            A = view(A, od1s...)
                        end
                        RasterPlot(), _prepare(_subsample(A, max_res))
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
    else
        first(st)
    end
end

# plot multispectral images as rgb
@userplot RGBPlot

@recipe function f(p::RGBPlot; bands=1:3, low=0.02, high=0.98)
    r = first(p.args)
    if !(typeof(r) <: AbstractRaster)
        error("Plotting RGB images is only implemented for AbstractRasters.")
    end
    if !hasdim(r, Band)
        error("Raster must have a band dimension.")
    end
    if !(0 <= low <= 1) || !(0 <= high <= 1)
        error("'low' and 'high' have to be between 0 and 1 (inclusive).")
    end
    img = Float32.(copy(r[Band([bands...])]))
    normalize!(img, low, high)
    img = permutedims(img, (Band, X, Y))
    img = DD.reorder(img, X=>DD.ForwardOrdered, Y=>DD.ForwardOrdered)
    plot_image = ImageCore.colorview(RGB, img)

    yguide, xguide = DD.label(dims(img, (Y,X)))

    y, x = map(Rasters._prepare, dims(img, (Y,X)))

    rdt = DD.refdims_title(img; issingle=true)
    :title --> (rdt === "" ? Rasters._maybename(img) : Rasters._maybename(img) * "\n" * rdt)
    :xguide --> xguide
    :xrotation --> 45
    :yguide --> yguide
    :tick_direction --> :out
    :framestyle --> :box

    if all(d -> lookup(d) isa Mapped, (x, y))
        :xlims --> mappedbounds(x)
        :ylims --> mappedbounds(y)
        :aspect_ratio --> :equal
    else
        :xlims --> bounds(img, x)
        :ylims --> bounds(img, y)
        bnds = bounds(img, (X, Y))
        # # TODO: Rethink this....
        s1, s2 = map(((l, u),) -> (u - l), bnds) ./ (size(img, X), size(img, Y))
        square_pixels = s2 / s1
        :aspect_ratio --> square_pixels
    end

    :seriestype --> :image
    :yflip --> false
    DD.index(x), DD.index(y), plot_image'
end


# Plots.jl heatmaps pixels are centered.
# So we should center the index, and use the projected value.
_prepare(d::Dimension) = d |> _maybe_shift |> _maybe_mapped
# Convert arrays to a consistent missing value and Forward array order
function _prepare(A::AbstractRaster)
    reorder(A, DD.ForwardOrdered) |>
    a -> permutedims(a, DD.commondims(>:, (ZDim, YDim, XDim, TimeDim, Dimension), dims(A)))# |>
end

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

_maybename(A) = _maybename(name(A))
_maybename(n::Name{N}) where N = _maybename(N)
_maybename(n::NoName) = ""
_maybename(n::Symbol) = string(n)
_maybename(n::AbstractString) = n

_maybe_replace_missing(A::AbstractArray{<:AbstractFloat}) = replace_missing(A, eltype(A)(NaN))
_maybe_replace_missing(A) = A

_maybe_shift(d) = _maybe_shift(sampling(d), d)
_maybe_shift(::Intervals, d) = DD.maybeshiftlocus(Center(), d)
_maybe_shift(sampling, d) = d

_maybe_mapped(dims::Tuple) = map(_maybe_mapped, dims)
_maybe_mapped(dim::Dimension) = _maybe_mapped(lookup(dim), dim)
_maybe_mapped(lookup::LookupArray, dim::Dimension) = dim
_maybe_mapped(lookup::Projected, dim::Dimension) = _maybe_mapped(mappedcrs(lookup), dim)
_maybe_mapped(::Nothing, dim::Dimension) = dim
_maybe_mapped(::GeoFormat, dim::Dimension) = convertlookup(Mapped, dim)

# We don't show the Band label for a single-band raster,
# it's just not interesting information.
function DD.refdims_title(refdim::Band; issingle=false)
    if issingle
        ""
    else
        string(name(refdim), ": ", DD.refdims_title(lookup(refdim), refdim))
    end
end

function eachband_view(r::Raster)
    nbands = size(r, Band)
    return (view(r, Band(n)) for n in 1:nbands)
end

function normalize!(raster, low=0.1, high=0.9)
    if !hasdim(raster, Band)
        l = quantile(skipmissing(raster), low)
        h = quantile(skipmissing(raster), high)
        raster .-= l
        raster ./= h - l + eps(float(eltype(raster)))
        raster .= clamp.(raster, zero(eltype(raster)), one(eltype(raster)))
    else
        for band in eachband_view(raster)
            l = quantile(skipmissing(band), low)
            h = quantile(skipmissing(band), high)
            band .-= l
            band ./= h - l + eps(float(eltype(raster)))
            band .= clamp.(band, zero(eltype(raster)), one(eltype(raster)))
        end
    end
    return raster
end