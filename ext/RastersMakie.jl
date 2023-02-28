module RastersMakie

# utility macro for backwards compatibility
macro _using(args...)
    @static if !isdefined(Base, :get_extension) # julia < 1.9
        Expr(:using, args)
    else # julia ≥ 1.9
        Expr(:using, Expr(:., :., :., args))
    end
end

@static if isdefined(Base, :get_extension) # julia < 1.9
    using Makie, Rasters
else    
    using ..Makie
    using ..Rasters
end

using Rasters.DimensionalData
using Rasters.MakieCore

# first, some Makie utils which require Makie types
MakieCore.plottype(::AbstractRaster{<: Makie.Colors.Colorant, 2}) = MakieCore.Image
function MakieCore.convert_arguments(::MakieCore.SurfaceLike, raw_raster::AbstractRaster{<: Makie.Colors.Colorant, 2})
    ds = DD._fwdorderdims(raw_raster)
    A = permutedims(raw_raster, ds)
    x, y = dims(A)
    xs, ys, zs = DD._withaxes(x, y, (A))
    return (xs, ys, collect(zs))
end

# now, the "full" plot-func

"""
    _balance_grid(nplots)

Returns a tuple `(ncols, nrows)`, defining a grid which can hold `nplots` items.
"""
function _balance_grid(nplots)
    ncols = (nplots - 1) ÷ ceil(Int, sqrt(nplots)) + 1
    nrows = (nplots - 1) ÷ ncols + 1
    return nrows, ncols
end


# The all-inclusive plotting function for a 2D raster

function Rasters.rplot(position::GridPosition, raster::AbstractRaster{T,2,<:Tuple{D1,D2}};
    plottype = Makie.Heatmap,
    axistype = Makie.Axis,
    X=X, Y=Y,
    draw_colorbar = true,
    colorbar_position = Makie.Right(),
    colorbar_padding = Makie.automatic,
    title = Makie.automatic,
    xlabel = Makie.automatic,
    ylabel = Makie.automatic,
    colorbarlabel = Makie.automatic,
    nan_color = (:brown, 0.02),
    colormap = nothing,
    colorrange = Makie.automatic,
    kw_attributes...
    ) where {T,D1<:Rasters.SpatialDim,D2<:Rasters.SpatialDim}

    # handle extra kwargs
    attributes = merge(
        Attributes(kw_attributes),
        Attributes(;
            Colorbar = (; label = colorbarlabel,),
            nan_color,
        )
    )

    isnothing(colormap) && (colormap = get(attributes, :colormap, :viridis))

    # x and y labels
    ylabel_str, xlabel_str = Rasters.label(DimensionalData.dims(raster))
    xlabel isa Makie.Automatic && (xlabel = xlabel_str)
    ylabel isa Makie.Automatic && (ylabel = ylabel_str)

    # colorbar label
    colorbarlabel isa Makie.Automatic && (colorbarlabel = ""#=string(DimensionalData.name(raster))=#)
    
    # title
    if title == Makie.automatic
        _rdt = DimensionalData.refdims_title(raster; issingle = true)
        title = (_rdt === "" ? Rasters._maybename(raster) : Rasters._maybename(raster) * "\n" * _rdt)
    end
    
    local axis, plot
    # actually plot
    with_theme(attributes) do
        axis = axistype(position; 
            title, xlabel, ylabel 
        )
        # plot to the axis with the specified plot type
        plot = plot!(plottype, axis, raster; colormap, colorrange, nan_color)

        if draw_colorbar

            colorbar = Colorbar(
                position.layout[position.span.rows, position.span.cols, colorbar_position],
                plot;
                label = colorbarlabel,
            )

            colorbar_padding = if colorbar_padding isa Makie.Automatic
                    if colorbar_position in (Makie.Top(), Makie.Bottom())
                        to_value(colorbar.layoutobservables.computedbbox[].widths[2])
                    else
                        to_value(colorbar.layoutobservables.computedbbox[].widths[1])
                    end
                else
                    to_value(colorbar_padding)
                end

            colorbar.alignmode[] = if colorbar_position == Makie.Top()
                Makie.Mixed(Makie.GridLayoutBase.RectSides{Union{Nothing, Float32, Makie.GridLayoutBase.Protrusion}}(nothing, nothing, colorbar_padding, 0f0))
            elseif colorbar_position == Makie.Left()
                Makie.Mixed(Makie.GridLayoutBase.RectSides{Union{Nothing, Float32, Makie.GridLayoutBase.Protrusion}}(0f0, colorbar_padding, nothing, nothing))
            elseif colorbar_position == Makie.Bottom()
                Makie.Mixed(Makie.GridLayoutBase.RectSides{Union{Nothing, Float32, Makie.GridLayoutBase.Protrusion}}(nothing, nothing, 0f0, colorbar_padding))
            elseif colorbar_position == Makie.Right()
                Makie.Mixed(Makie.GridLayoutBase.RectSides{Union{Nothing, Float32, Makie.GridLayoutBase.Protrusion}}(colorbar_padding, 0f0, nothing, nothing))
            else
                @error "The colorbar position `$(colorbar_position)` was not recognized.  Please pass one of `Makie.Top(), Makie.Bottom(), Makie.Right(), Makie.Left()`."
            end
        end
    end
    return Makie.AxisPlot(axis, plot)
end

function Rasters.rplot(raster::AbstractRaster{T, 2}; kwargs...) where T
    figure = isempty(kwargs) ? Figure() : with_theme(Figure, Attributes(kwargs))
    axis, plot = Rasters.rplot(figure[1, 1], raster; kwargs...)
    return Makie.FigureAxisPlot(figure, axis, plot)
end

function Rasters.rplot(gp::GridPosition, raster::AbstractRaster{T, 3}; ncols = Makie.automatic, nrows = Makie.automatic, kwargs...) where T

    nrows, ncols = if ncols isa Makie.Automatic && nrows isa Makie.Automatic
        _balance_grid(size(raster, 3))
    elseif ncols isa Int && nrows isa Int
        @assert ncols * nrows ≥ size(raster, 3)
        nrows, ncols
    else
        @error("The provided combination of `ncols::$(typeof(ncols)) and nrows::$(typeof(nrows)) is unsupported.  Please either set both to `Makie.automatic` or provide integer values.")
    end

    layout = GridLayout(gp, nrows, ncols)

    for (i, band) in enumerate(axes(raster, 3))
        ax, plt = Rasters.rplot(layout[fldmod1(i, ncols)...], view(raster, :, :, band); kwargs...)
        if fldmod1(i, ncols)[2] != 1
            hideydecorations!(ax, label = true, ticklabels = true, ticks = false, grid = false, minorgrid = false, minorticks = false)
        end
        if fldmod1(i, ncols)[1] != nrows
            hidexdecorations!(ax, label = true, ticklabels = true, ticks = false, grid = false, minorgrid = false, minorticks = false)
        end
    end

    return layout
end

function Rasters.rplot(gp::GridPosition, stack::RasterStack; ncols = Makie.automatic, nrows = Makie.automatic, colormap = nothing, colorrange = Makie.automatic, link_colorrange = false, link_axes = true, kwargs...) where T

    @assert (length(size(stack)) == 2 || size(stack, 3) == 1)

    nrows, ncols = if ncols isa Makie.Automatic && nrows isa Makie.Automatic
        _balance_grid(length(propertynames(stack)))
    elseif ncols isa Int && nrows isa Int
        @assert ncols * nrows ≥ length(propertynames(stack))
        nrows, ncols
    else
        @error("The provided combination of `ncols::$(typeof(ncols)) and nrows::$(typeof(nrows)) is unsupported.  Please either set both to `Makie.automatic` or provide integer values.")
    end

    layout = GridLayout(gp, nrows, ncols)

    axs = [] # avoid ambiguity with Base.axes
    plots = []

    Makie.broadcast_foreach(collect(enumerate(propertynames(stack))), colormap, colorrange) do (i, band), cmap, crange
        raster = getproperty(stack, band)
        if length(size(raster)) == 2
            raster = raster
        elseif length(size(raster)) == 3
            if size(raster, 3) == 1
                raster = view(raster, :, :, 1)
            else
                @error "You can't plot a RasterStack of 3-D rasters using `rplot`.  Please provide a stack of 2D rasters instead, or 3D rasters with a singleton third dimension."
            end
        else
            @error "`rplot` cannot plot a Raster of dimension $(size(raster)).  Please provide a stack of 2D rasters instead."
        end

        ax, plt = Rasters.rplot(layout[fldmod1(i, ncols)...], raster; colormap = cmap, colorrange = crange, kwargs...)

        if fldmod1(i, ncols)[2] != 1
            hideydecorations!(ax, label = true, ticklabels = true, ticks = false, grid = false, minorgrid = false, minorticks = false)
        end
        if fldmod1(i, ncols)[1] != nrows
            hidexdecorations!(ax, label = true, ticklabels = true, ticks = false, grid = false, minorgrid = false, minorticks = false)
        end

        push!(axs, ax)
        push!(plots, plt)
    end

    if link_axes
        linkaxes!(axs...)
    end

    if link_colorrange
        lift(getproperty.(plots, :colorrange)...; ignore_equal_values = true) do cranges...
            new_colorrange = (minimum(first.(cranges)), maximum(last.(cranges)))
            setproperty!.(plots, :colorrange, (new_colorrange,))
        end
    end

    return layout
end

function Rasters.rplot(raster::AbstractRaster{T, 3}; colormap = nothing, colorrange = Makie.automatic, kwargs...) where T
    figure = isempty(kwargs) ? Figure() : with_theme(Figure, Attributes(kwargs))
    layout = Rasters.rplot(figure[1, 1], raster; colormap, colorrange, kwargs...)
    # if draw_title
    #     Label(layout[0, 1:Makie.ncols(layout)], raster_title; fontsize = get(figure.scene.attributes, (:Axis :titlesize), 16), font = get(figure.scene.attributes, (:Axis, :titlefont), :bold))
    # end
    return figure
end


function Rasters.rplot(raster::RasterStack; kwargs...)
    figure = isempty(kwargs) ? Figure() : with_theme(Figure, Attributes(kwargs))
    layout = Rasters.rplot(figure[1, 1], raster; kwargs...)
    # if draw_title
    #     Label(layout[0, 1:Makie.ncols(layout)], raster_title; fontsize = get(figure.scene.attributes, (:Axis :titlesize), 16), font = get(figure.scene.attributes, (:Axis, :titlefont), :bold))
    # end
    return figure
end

end

# fig, ax1, plt1 = rplot(clamp.(raster_2012, 0, 20))
# ax2, plt2 = rplot(fig[1, 2], clamp.(raster_2012, 0, 20); draw_colorbar = true, Axis = (; aspect = DataAspect()))
# fig
