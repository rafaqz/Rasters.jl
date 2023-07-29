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

# For now, these are just basic conversions from 2D rasters to surface-like (surface, heatmap) plot types.
# Once Makie supports figure level recipes, we should integrate those.

# First, define default plot types for Rasters
MakieCore.plottype(raster::AbstractRaster{<:Union{Missing,Real},1}) = MakieCore.Lines
MakieCore.plottype(raster::AbstractRaster{<:Union{Missing,Real},2}) = MakieCore.Heatmap
# 3d rasters are a little more complicated - if dim3 is a singleton, then heatmap, otherwise volume
function MakieCore.plottype(raster::AbstractRaster{<:Union{Missing,Real},3})
    if size(raster, 3) == 1
        MakieCore.Heatmap
    else
        MakieCore.Volume
    end
end
# plot rasters of ColorTypes as images
# define the correct plottype
MakieCore.plottype(::AbstractRaster{<:ColorTypes.Colorant,2}) = MakieCore.Image

# then, define how they are to be converted to plottable data
function MakieCore.convert_arguments(PB::MakieCore.PointBased, raster::AbstractRaster{<:Union{Real,Missing},1})
    A = _prepare_makie(raster, ds)
    return MakieCore.convert_arguments(PB, A, index(ds))
end
# allow 3d rasters to be plotted as volumes
function MakieCore.convert_arguments(
    ::MakieCore.VolumeLike, raster::AbstractRaster{<:Union{Real,Missing},3,<:Tuple{D1,D2,D3}}
) where {D1<:SpatialDim,D2<:SpatialDim,D3<:SpatialDim}
    A = _prepare_makie(raster)
    xs, ys, zs = lookup(A)
    return (xs, ys, zs, parent(A))
end
# allow plotting 3d rasters with singleton third dimension (basically 2d rasters)
function MakieCore.convert_arguments(x::MakieCore.ConversionTrait, raster::AbstractRaster{<:Union{Real,Missing},3})
    D = _series_dim(raster)
    nplots = size(raster, D)
    if nplots > 1
        # Plot as a RasterSeries
        return MakieCore.convert_arguments(x, slice(raster, D))
    else
        return MakieCore.convert_arguments(x, view(raster, rebuild(D, 1)))
    end
end
function MakieCore.convert_arguments(
    ::MakieCore.SurfaceLike, raster::AbstractRaster{<:Union{Missing,Real},2,<:Tuple{D1,D2}}
) where {D1<:SpatialDim,D2<:SpatialDim}
    A = _prepare_makie(raster)
    xs, ys = lookup(A)
    return (xs, ys, parent(A))
end
function MakieCore.convert_arguments(
    ::MakieCore.SurfaceLike, raster::AbstractRaster{<:ColorTypes.Colorant,2,<:Tuple{D1,D2}}
) where {D1<:SpatialDim,D2<:SpatialDim}
    A = _prepare_makie(raster)
    x, y = dims(A)
    xs, ys, vs = DD._withaxes(x, y, A)
    return (xs, ys, parent(vs))
end
function MakieCore.convert_arguments(
    ::MakieCore.DiscreteSurface, raster::AbstractRaster{<:Union{Missing,Real},2,<:Tuple{D1,D2}}
) where {D1<:SpatialDim,D2<:SpatialDim}
    A = _prepare_makie(raster)
    xs, ys = map(_lookup_edges, lookup(A))
    return (xs, ys, parent(A))
end
function MakieCore.convert_arguments(
    ::MakieCore.DiscreteSurface,raster::AbstractRaster{<:ColorTypes.Colorant,2,<:Tuple{D1,D2}}
) where {D1<:SpatialDim,D2<:SpatialDim}
    A = _prepare_makie(raster)
    xs, ys = map(_lookup_edges, lookup(A, (X, Y)))
    return (xs, ys, parent(A))
end
function MakieCore.convert_arguments(x::MakieCore.ConversionTrait, series::AbstractRasterSeries)
    return MakieCore.convert_arguments(x, first(series))
end
# fallbacks with descriptive error messages
MakieCore.convert_arguments(t::MakieCore.ConversionTrait, r::AbstractRaster) =
    _makie_not_implemented_error(t, r)

function _makie_not_implemented_error(t, r::AbstractRaster{T,N}) where {T,N}
    @error """
    We don't currently support plotting Rasters of eltype $T, with $N dimensions, as
    a $t type plot in Makie.jl.

    In order to plot, please provide a 2-dimensional slice of your raster contining 
    only `Real` and `Missing` values.

    For example, in a 3-dimensional raster, this would look like:

    ```julia
    myraster = Raster(...)
    heatmap(myraster)        # errors
    heatmap(myraster[Ti(1)]) # use some index to subset, this works!
    ```
    """
end

_prepare_makie(A) = 
    _reorder(read(_missing_or_float32.(replace_missing(A; missingval=NaN32))))

# initial definitions of `rplot`, to get around the extension package availability question

function rplot(args...)
    @error("Please load `Makie.jl` and then call this function. If Makie is loaded, then you can't call `rplot` with no arguments!")
end

# define the theme

# this function is defined so that we can override style_rasters in RastersMakieExt
function style_rasters end

function color_rasters()
    return MakieCore.Attributes(
        colormap = :plasma,
    )
end

function theme_rasters()
    return merge(style_rasters(), color_rasters())
end


##################################################################################
# Utils

_missing_or_float32(num::Number) = Float32(num)
_missing_or_float32(::Missing) = missing

function _lookup_edges(l::LookupArray)
    l = if l isa AbstractSampled 
        set(l, Intervals())
    else
        set(l, Sampled(; sampling=Intervals()))
    end
    if l == 1
        return [bounds(l)...]
    else
        ib = intervalbounds(l)
        if order(l) isa ForwardOrdered
            edges = first.(ib)
            push!(edges, last(last(ib)))
        else
            edges = last.(ib)
            push!(edges, first(last(ib)))
        end
        return edges
    end
end

# Plots.jl heatmaps pixels are centered.
# So we should center the index, and use the projected value.
_prepare(d::Dimension) = d |> _maybe_shift |> _maybe_mapped
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

function _series_dim(A)
    spatialdims = (X(), Y(), Z())
    last((dims(A, spatialdims)..., otherdims(A, spatialdims)...))
end
