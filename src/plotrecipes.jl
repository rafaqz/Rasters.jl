# Method specialisation singletons.
struct RasterPlot end
struct RasterZPlot end
# We only look at arrays with X, Y, Z dims here.
# Otherwise they fall back to DimensionalData.jl recipes
@recipe function f(A::AbstractRaster)
    ddplot(A) = DimArray(A; dims=_maybe_mapped(dims(A)))
    # Resample AffineProjected
    max_res = get(plotattributes, :max_res, 1000)
    A = maybe_resample(A)
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

function maybe_resample(A)
    if all(hasdim(A, (X, Y))) 
        maybe_resample(first(lookup(A, (X, Y))), A)
    else
        return A
    end
end
maybe_resample(lookup, A) = A



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
    nplots = length(A)
    if nplots > 16 
        plotinds = round.(Int, 1:nplots//16:nplots)
        @info "Too many raster heatmaps: plotting 16 slices from $nplots"
        A = @views A[plotinds]
        nplots = length(plotinds)
    end
    # link --> :both
    # :colorbar := false
    titles = string.(index(A, dims(A, 1)))
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
                    A[i]
                else
                    :framestyle := :none
                    :legend := :none
                    []
                end
            end
        end
    end
end

# Makie.jl recipes

# For now, these are just basic conversions from 2D rasters to surface-like (surface, heatmap) plot types.
# Once Makie supports figure level recipes, we should integrate those.
    
# First, define default plot types for Rasters
MakieCore.plottype(raw_raster::AbstractRaster{<: Union{Missing, Real}, 1}) = MakieCore.Lines
MakieCore.plottype(raw_raster::AbstractRaster{<: Union{Missing, Real}, 2}) = MakieCore.Heatmap
# 3d rasters are a little more complicated - if dim3 is a singleton, then heatmap, otherwise volume
function MakieCore.plottype(raw_raster::AbstractRaster{<: Union{Missing, Real}, 3})
    if size(raw_raster, 3) == 1
        MakieCore.Heatmap
    else
        MakieCore.Volume
    end
end

missing_or_float32(num::Number) = Float32(num)
missing_or_float32(::Missing) = missing

# then, define how they are to be converted to plottable data
function MakieCore.convert_arguments(PB::MakieCore.PointBased, raw_raster::AbstractRaster{<: Union{Missing, Real}, 1})
    z = map(Rasters._prepare, dims(raw_raster))
    return MakieCore.convert_arguments(PB, parent(Float32.(replace_missing(missing_or_float32.(raw_raster), missingval = NaN32))), index(z))
end

function MakieCore.convert_arguments(::MakieCore.SurfaceLike, raw_raster::AbstractRaster{<: Union{Missing, Real}, 2})
    raster = replace_missing(missing_or_float32.(raw_raster), missingval = NaN32)
    ds = DD._fwdorderdims(raster)
    A = permutedims(raster, ds)
    x, y = dims(A)
    xs, ys, vs = DD._withaxes(x, y, (A))
    return (xs, ys, collect(vs))
end

function __edges(v::AbstractVector)
    l = length(v)
    if l == 1
        return [v[begin] - 0.5, v[begin] + 0.5]
    else
        # Equivalent to
        # mids = 0.5 .* (v[1:end-1] .+ v[2:end])
        # borders = [2v[1] - mids[1]; mids; 2v[end] - mids[end]]
        borders = [0.5 * (v[max(begin, i)] + v[min(end, i+1)]) for i in (firstindex(v) - 1):lastindex(v)]
        borders[1] = 2borders[1] - borders[2]
        borders[end] = 2borders[end] - borders[end-1]
        return borders
    end
end

function MakieCore.convert_arguments(::MakieCore.DiscreteSurface, raw_raster::AbstractRaster{<: Union{Missing, Real}, 2})
    raster = replace_missing(missing_or_float32.(raw_raster), missingval = NaN32)
    ds = DD._fwdorderdims(raster)
    A = permutedims(raster, ds)
    x, y = dims(A)
    xs, ys, vs = DD._withaxes(x, y, (A))
    return (__edges(xs), __edges(ys), collect(vs))
end

# allow plotting 3d rasters with singleton third dimension (basically 2d rasters)
function MakieCore.convert_arguments(sl::MakieCore.SurfaceLike, raw_raster_with_missings::AbstractRaster{<: Union{Real, Missing}, 3})
    @assert size(raw_raster_with_missings, 3) == 1
    return MakieCore.convert_arguments(sl, raw_raster_with_missings[:, :, begin])
end

# allow 3d rasters to be plotted as volumes
function MakieCore.convert_arguments(::MakieCore.VolumeLike, raw_raster_with_missings::AbstractRaster{<: Union{Real, Missing}, 3})
    raster = replace_missing(missing_or_float32.(raw_raster), missingval = NaN32)
    ds = DD._fwdorderdims(raster)
    A = permutedims(raster, ds)
    x, y, z = dims(A)
    xs, ys, zs, vs = DD._withaxes(x, y, z, A)
    return (xs, ys, zs, collect(vs))
end

# plot rasters of ColorTypes as images
# define the correct plottype
MakieCore.plottype(::AbstractRaster{<: ColorTypes.Colorant, 2}) = MakieCore.Image

function MakieCore.convert_arguments(::MakieCore.SurfaceLike, raw_raster::AbstractRaster{<: ColorTypes.Colorant, 2})
    ds = DD._fwdorderdims(raw_raster)
    A = permutedims(raw_raster, ds)
    x, y = dims(A)
    xs, ys, vs = DD._withaxes(x, y, (A))
    return (xs, ys, collect(vs))
end

function MakieCore.convert_arguments(::MakieCore.DiscreteSurface, raw_raster::AbstractRaster{<: ColorTypes.Colorant, 2})
    ds = DD._fwdorderdims(raw_raster)
    A = permutedims(raw_raster, ds)
    x, y = dims(A)
    xs, ys, vs = DD._withaxes(x, y, (A))
    return (__edges(xs), __edges(ys), collect(vs))
end
            
# fallbacks with descriptive error messages
MakieCore.convert_arguments(::MakieCore.SurfaceLike, ::AbstractRaster{<: Real, Dim}) = @error """
            We don't currently support plotting Rasters of dimension $Dim in Makie.jl. in surface/heatmap-like plots.
            
            In order to plot, please provide a 2-dimensional slice of your raster.
            For example, in a 3-dimensional raster, this would look like:
            ```julia
            myraster = Raster(...)
            heatmap(myraster)          # errors
            heatmap(myraster[:, :, 1]) # use some index to subset, this works!
            ```
            """
            

# initial definitions of `rplot`, to get around the extension package availability question

function rplot() 
    @error("Please load `Makie.jl` and then call this function.  If Makie is loaded, then you can't call `rplot` with no arguments!")
end

# define the theme

# this function is defined so that we can override style_rasters in RastersMakie.jl
function __style_rasters()
    return MakieCore.Attributes(
        Axis = (
            xtickalign = 1.0,
            ytickalign = 1.0,
            xticklabelrotation = -π/4,
            xticklabelsize = 14,
            yticklabelsize = 14,
            # aspect = DataAspect(),
        ),

        Colorbar = (
            ticklabelsize = 11,
            tickalign = 1.0,
        ),
    )
end

function style_rasters end # defined in ../ext/RastersMakie

function color_rasters()
    return MakieCore.Attributes(
        colormap = :plasma,
    )
end

function theme_rasters()
    return merge(style_rasters(), color_rasters())
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

function _balance_grid(nplots)
    ncols = (nplots - 1) ÷ ceil(Int, sqrt(nplots)) + 1
    nrows = (nplots - 1) ÷ ncols + 1
    return ncols, nrows
end
