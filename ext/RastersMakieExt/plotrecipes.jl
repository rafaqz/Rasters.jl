function Rasters.style_rasters()
    return merge(
        Rasters.__style_rasters(),
        Attributes(
            Axis = (; aspect = DataAspect()),
        )
    )
end

function lift_layer(r::Observable, inds...)
    return lift(lift_layer, r, inds...)
end
lift_layer(r::Raster, inds...) = getindex(r, inds...)
lift_layer(rs::RasterStack, ind::Symbol) = getproperty(rs, ind)

# The all-inclusive plotting function for a 2D raster
"""
    Rasters.rplot([position::GridPosition], raster; kw_args...)

`raster` may be a `Raster` (of 2 or 3 dimensions) or a `RasterStack` whose underlying rasters are 2 dimensional, or 3-dimensional with a singleton (length-1) third dimension.

## Keyword arguments
- `plottype = Makie.Heatmap`: The type of plot.  Can be any Makie plot type which accepts a `Raster`; in practice, `Heatmap`, `Contour`, `Contourf` and `Surface` are the best bets.
- `axistype = Makie.Axis`: The type of axis.  This can be an `Axis`, `Axis3`, `LScene`, or even a `GeoAxis` from GeoMakie.jl.
- `X=X`: The X dimension of the raster.
- `Y=Y`: The Y dimension of the raster.
- `draw_colorbar = true`: Whether to draw a colorbar for the axis or not.
- `colorbar_position = Makie.Right()`: Indicates which side of the axis the colorbar should be placed on.  Can be `Makie.Top()`, `Makie.Bottom()`, `Makie.Left()`, or `Makie.Right()`.
- `colorbar_padding = Makie.automatic`: The amound of padding between the colorbar and its axis.  If `automatic`, then this is set to the width of the colorbar.
- `title = Makie.automatic`: The titles of each plot.  If `automatic`, these are set to the name of the band.
- `xlabel = Makie.automatic`: The x-label for the axis.  If `automatic`, set to the dimension name of the X-dimension of the raster.
- `ylabel = Makie.automatic`: The y-label for the axis.  If `automatic`, set to the dimension name of the Y-dimension of the raster.
- `colorbarlabel = ""`: Usually nothing, but here if you need it.  Sets the label on the colorbar.
- `colormap = nothing`: The colormap for the heatmap.  This can be set to a vector of colormaps (symbols, strings, `cgrad`s) if plotting a 3D raster or RasterStack.
- `colorrange = Makie.automatic`: The colormap for the heatmap.  This can be set to a vector of `(low, high)` if plotting a 3D raster or RasterStack.
- `nan_color = :transparent`: The color which `NaN` values should take.  Default to transparent.
"""
function Rasters.rplot(position::GridPosition, raster::Union{AbstractRaster{T,2,<:Tuple{D1,D2}}, Observable{<: AbstractRaster{T,2,<:Tuple{D1,D2}}}};
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
    nan_color = :transparent,
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

    val_raster = Makie.to_value(raster)

    isnothing(colormap) && (colormap = get(attributes, :colormap, :viridis))

    # x and y labels
    xlabel_str, ylabel_str = Rasters.label(DimensionalData.dims(val_raster))
    xlabel isa Makie.Automatic && (xlabel = xlabel_str)
    ylabel isa Makie.Automatic && (ylabel = ylabel_str)

    # colorbar label
    colorbarlabel isa Makie.Automatic && (colorbarlabel = ""#=string(DimensionalData.name(raster))=#)
    
    # title
    if title == Makie.automatic
        _rdt = DimensionalData.refdims_title(val_raster; issingle = true)
        title = (_rdt === "" ? Rasters._maybename(val_raster) : Rasters._maybename(val_raster) * "\n" * _rdt)
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
                Makie.Mixed(bottom = colorbar_padding, #= top = 0 =#)
            elseif colorbar_position == Makie.Left()
                Makie.Mixed(right = colorbar_padding, #= left = 0 =#)
            elseif colorbar_position == Makie.Bottom()
                Makie.Mixed(top = colorbar_padding, #= bottom = 0 =#)
            elseif colorbar_position == Makie.Right()
                Makie.Mixed(left = colorbar_padding, #= right = 0 =#)
            else
                @error "The colorbar position `$(colorbar_position)` was not recognized.  Please pass one of `Makie.Top(), Makie.Bottom(), Makie.Right(), Makie.Left()`."
            end
        end
    end
    return Makie.AxisPlot(axis, plot)
end

function Rasters.rplot(gp::GridPosition, raster::Union{AbstractRaster{T, 3}, Observable{<: AbstractRaster{T, 3}}}; ncols = Makie.automatic, nrows = Makie.automatic, kwargs...) where T

    val_raster = Makie.to_value(raster)

    nrows, ncols = if ncols isa Makie.Automatic && nrows isa Makie.Automatic
        Rasters._balance_grid(size(val_raster, 3))
    elseif ncols isa Int && nrows isa Int
        @assert ncols * nrows ≥ size(val_raster, 3)
        nrows, ncols
    else
        @error("The provided combination of `ncols::$(typeof(ncols)) and nrows::$(typeof(nrows)) is unsupported.  Please either set both to `Makie.automatic` or provide integer values.")
    end

    layout = GridLayout(gp, nrows, ncols)

    for (i, layer) in enumerate(axes(val_raster, 3))
        ax, plt = Rasters.rplot(layout[fldmod1(i, ncols)...], lift_layer(raster, :, :, layer); kwargs...)
        if fldmod1(i, ncols)[2] != 1
            hideydecorations!(ax, label = true, ticklabels = true, ticks = false, grid = false, minorgrid = false, minorticks = false)
        end
        if fldmod1(i, ncols)[1] != nrows
            hidexdecorations!(ax, label = true, ticklabels = true, ticks = false, grid = false, minorgrid = false, minorticks = false)
        end
    end

    return layout
end

function Rasters.rplot(gp::GridPosition, stack::Union{RasterStack, Observable{<: RasterStack}}; ncols = Makie.automatic, nrows = Makie.automatic, colormap = nothing, colorrange = Makie.automatic, link_colorrange = false, link_axes = true, axis = (;), kwargs...)

    val_stack = Makie.to_value(stack)

    @assert (length(size(val_stack)) == 2 || size(val_stack, 3) == 1)

    nrows, ncols = if ncols isa Makie.Automatic && nrows isa Makie.Automatic
        Rasters._balance_grid(length(propertynames(val_stack)))
    elseif ncols isa Int && nrows isa Int
        @assert ncols * nrows ≥ length(propertynames(val_stack))
        nrows, ncols
    else
        @error("The provided combination of `ncols::$(typeof(ncols)) and nrows::$(typeof(nrows)) is unsupported.  Please either set both to `Makie.automatic` or provide integer values.")
    end

    layout = GridLayout(gp, nrows, ncols)

    axs = [] # avoid ambiguity with Base.axes
    plots = []

    Makie.broadcast_foreach(collect(enumerate(propertynames(val_stack))), colormap, colorrange) do (i, layer), cmap, crange
        raster = lift_layer(stack, layer)
        if length(size(to_value(raster))) == 2
            raster = raster
        elseif length(size(to_value(raster))) == 3
            if size(to_value(raster), 3) == 1
                raster = lift_layer(raster, Rasters.Band(1))
            else
                @error "You can't plot a RasterStack of 3-D rasters using `rplot`.  Please provide a stack of 2D rasters instead, or 3D rasters with a singleton third dimension."
            end
        else
            @error "`rplot` cannot plot a Raster of dimension $(size(raster)).  Please provide a stack of 2D rasters instead."
        end

        ax, plt = Rasters.rplot(layout[fldmod1(i, ncols)...], raster; colormap = cmap, colorrange = crange, axis = axis, kwargs...)

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

function Rasters.rplot(raster::Union{AbstractRaster{T, 2}, Observable{<: AbstractRaster{T, 2}}}; figure = (;), kwargs...) where T
    figure = Figure(; figure...)
    axis, plot = Rasters.rplot(figure[1, 1], raster; kwargs...)
    return Makie.FigureAxisPlot(figure, axis, plot)
end

function Rasters.rplot(raster::Union{AbstractRaster{T, 3}, Observable{<: AbstractRaster{T, 3}}}; figure = (;), colormap = nothing, colorrange = Makie.automatic, kwargs...) where T
    figure = Figure(; figure...)
    layout = Rasters.rplot(figure[1, 1], raster; colormap, colorrange, kwargs...)
    # if draw_title
    #     Label(layout[0, 1:Makie.ncols(layout)], raster_title; fontsize = get(figure.scene.attributes, (:Axis :titlesize), 16), font = get(figure.scene.attributes, (:Axis, :titlefont), :bold))
    # end
    return figure
end

function Rasters.rplot(raster::Union{RasterStack, Observable{<: RasterStack}}; figure = (;), colormap = nothing, colorrange = Makie.automatic, kwargs...)
    figure = Figure(; figure...)
    layout = Rasters.rplot(figure[1, 1], raster; colormap, colorrange, kwargs...)
    # if draw_title
    #     Label(layout[0, 1:Makie.ncols(layout)], raster_title; fontsize = get(figure.scene.attributes, (:Axis :titlesize), 16), font = get(figure.scene.attributes, (:Axis, :titlefont), :bold))
    # end
    return figure
end
