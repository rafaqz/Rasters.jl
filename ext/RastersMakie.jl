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
    xlabel isa Makie.Automatic && (xlabel = xlabel_str)
    ylabel isa Makie.Automatic && (ylabel = ylabel_str)

    # colorbar label
    colorbarlabel isa Makie.Automatic && (colorbarlabel = string(DimensionalData.name(raster)))
    
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
        plot = plot!(plottype, axis, raster; nan_color)

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
                Makie.Outside(0, 0, colorbar_padding, 0)
            elseif colorbar_position == Makie.Left()
                Makie.Outside(0, colorbar_padding, 0, 0)
            elseif colorbar_position == Makie.Bottom()
                Makie.Outside(0, 0, 0, colorbar_padding)
            elseif colorbar_position == Makie.Right()
                Makie.Outside(colorbar_padding, 0, 0, 0)
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
        Rasters.rplot(layout[fldmod1(i, ncols)...], view(raster, :, :, band); kwargs...)
    end

    return layout
end

function Rasters.rplot(raster::AbstractRaster{T, 3}; kwargs...) where T
    figure = isempty(kwargs) ? Figure() : with_theme(Figure, Attributes(kwargs))
    layout = Rasters.rplot(figure[1, 1], raster; kwargs...)
    # if draw_title
    #     Label(layout[0, 1:Makie.ncols(layout)], raster_title; fontsize = get(figure.scene.attributes, (:Axis :titlesize), 16), font = get(figure.scene.attributes, (:Axis, :titlefont), :bold))
    # end
    return figure
end

end