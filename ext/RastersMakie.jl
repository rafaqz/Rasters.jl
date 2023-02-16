module RastersMakie

# utility macro for backwards compatibility
macro _using(args...)
    @static if !isdefined(Base, :get_extension) # julia < 1.9
        Expr(:using, args)
    else # julia ≥ 1.9
        Expr(:using, Expr(:., :., :., args))
    end
end

@_using Makie
@_using Rasters

using Rasters.DimensionalData

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
    return ncols, nrows
end


# The all-inclusive plotting function for a 2D raster

function rplot(position::GridPosition, raster::AbstractRaster{T,2,<:Tuple{D1,D2}};
    plottype = Makie.Heatmap,
    axistype = Makie.Axis,
    X=X, Y=Y,
    draw_colorbar = false,
    colorbar_position = Makie.Right(),
    colorbar_padding = 15,
    title = Makie.automatic,
    xlabel = Makie.automatic,
    ylabel = Makie.automatic,
    colorbarlabel = Makie.automatic,
    nan_color = (:brown, 0.02),
    replace_missing = false,
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

    # x and y labels
    ylabel_str, xlabel_str = Rasters.label(DimensionalData.dims(raster))
    xlabel == Makie.automatic && (xlabel = xlabel_str)
    ylabel == Makie.automatic && (ylabel = ylabel_str)

    # colorbar label
    colorbarlabel == Makie.automatic && (colorbarlabel = string(DimensionalData.name(raster)))
    
    # title
    if title == Makie.automatic
        _rdt = DimensionalData.refdims_title(raster; issingle = true)
        title = (_rdt === "" ? Rasters._maybename(raster) : Rasters._maybename(raster) * "\n" * _rdt)
    end
    
    local axis, plot
    # actually plot
    with_theme(attributes) do
        layout = GridLayout(position)
        axis = axistype(layout[1, 1]; 
            title, xlabel, ylabel 
        )
        # plot to the axis with the specified plot type
        plot = plot!(plottype, axis, raster; nan_color)

        if draw_colorbar
            current_colorbar_theme = get(Makie.current_default_theme(), :Colorbar, nothing)
            correct_layout_cell = if colorbar_position == Makie.Top()
                layout[0, 1]
            elseif colorbar_position == Makie.Left()
                layout[1, 0]
            elseif colorbar_position == Makie.Bottom()
                layout[2, 1]
            elseif colorbar_position == Makie.Right()
                layout[1, 2]
            else
                @error "The colorbar position `$(colorbar_position)` was not recognized.  Please pass one of `Makie.Top(), Makie.Bottom(), Makie.Right(), Makie.Left()`."
            end
                colorbar = Colorbar(
                #=position.layout[position.span.rows, position.span.cols, colorbar_position]=#
                correct_layout_cell, plot;
                label = colorbarlabel,
            )
        end
    end
    return Makie.AxisPlot(axis, plot)
end

function rplot(raster::AbstractRaster{T,2}; kwargs...) where T
    figure = isempty(kwargs) ? Figure() : with_theme(Figure, kwargs)
    axis, plot = rplot(figure[1, 1], raster; kwargs...)
    return Makie.FigureAxisPlot(figure, axis, plot)
end

end