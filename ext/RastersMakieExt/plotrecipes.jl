
const HIDE_DEC = (; label=true, grid=false, minorgrid=false, minorticks=false)

const SurfaceLikeCompat = isdefined(Makie, :SurfaceLike) ? Makie.SurfaceLike : Union{Makie.VertexGrid,Makie.CellGrid,Makie.ImageLike}

function Rasters.style_rasters()
    Makie.Attributes(
        Axis=(
            xtickalign=1.0,
            ytickalign=1.0,
            xticklabelrotation=-π/4,
            xticklabelsize=14,
            yticklabelsize=14,
            aspect=DataAspect(),
        ),
        Colorbar=(
            ticklabelsize=11,
            tickalign=1.0,
        ),
    )
end

function Rasters.color_rasters()
    return Makie.Attributes(
        colormap = :batlow,
    )
end

function lift_layer(r::Observable, inds...)
    return lift(lift_layer, r, inds...)
end
lift_layer(r::Raster, inds...) = getindex(r, inds...)
lift_layer(rs::RasterStack, ind::Symbol) = getproperty(rs, ind)
lift_layer(s::RasterSeries, inds...) = getindex(s, inds...)

# The all-inclusive plotting function for a 2D raster
function Rasters.rplot(position::GridPosition, raster::Union{AbstractRaster{T,2,<:Tuple{D1,D2}},Observable{<:AbstractRaster{T,2,<:Tuple{D1,D2}}}};
    axistype=Makie.Axis,
    X=XDim, Y=YDim, Z=ZDim,
    draw_colorbar=true,
    colorbar_position=Makie.Right(),
    colorbar_padding=Makie.automatic,
    title=Makie.automatic,
    xlabel=Makie.automatic,
    ylabel=Makie.automatic,
    colorbarlabel=Makie.automatic,
    nan_color=:transparent,
    colormap=nothing,
    colorrange=Makie.automatic,
    kw_attributes...
) where {T,D1<:Rasters.SpatialDim,D2<:Rasters.SpatialDim}

    val_raster = Makie.to_value(raster)

    # x and y labels
    xlabel_str, ylabel_str = Rasters.label(DimensionalData.dims(val_raster))
    if xlabel isa Makie.Automatic
        xlabel = xlabel_str
    end
    if ylabel isa Makie.Automatic
        ylabel = ylabel_str
    end

    # colorbar label
    if colorbarlabel isa Makie.Automatic
        colorbarlabel = "" # string(DimensionalData.name(raster))
    end

    # handle extra kw
    attributes = merge(
        Attributes(kw_attributes),
        Attributes(;
            Colorbar=(; label=colorbarlabel,),
            nan_color,
        )
    )

    if isnothing(colormap)
        colormap = get(attributes, :colormap, :batlow)
    end

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
        plot = heatmap!(axis, raster; colormap, colorrange, nan_color)

        if draw_colorbar
            layout = position.layout[position.span.rows, position.span.cols, colorbar_position]
            colorbar = Colorbar(layout, plot; label=colorbarlabel)

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
                @error "The colorbar position `$(colorbar_position)` was not recognized. Please pass one of `Makie.Top(), Makie.Bottom(), Makie.Right(), Makie.Left()`."
            end
        end
    end

    return Makie.AxisPlot(axis, plot)
end
function Rasters.rplot(
    gp::GridPosition, raster::Union{AbstractRaster{T,3},Observable{<:AbstractRaster{T,3}}};
    X=XDim, Y=YDim, Z=ZDim,
    kw...
) where T
    val_raster = Makie.to_value(raster)
    spatialdims = (XDim, YDim, ZDim)
    ds = (dims(val_raster, spatialdims)..., otherdims(val_raster, spatialdims)...)
    return Rasters.rplot(gp, slice(raster, last(ds)); X, Y, Z, kw...)
end
function Rasters.rplot(
    gp::GridPosition, 
    series::Union{<:AbstractRasterSeries{<:Any,1},Observable{<:AbstractRasterSeries{<:Any,1}}};
    ncols=Makie.automatic,
    nrows=Makie.automatic,
    kw...
)
    val_series = Makie.to_value(series)

    nrows, ncols, plotinds = if ncols isa Makie.Automatic && nrows isa Makie.Automatic
        _, plotinds, nplots = Rasters._maybe_thin_plots(val_series)
        nrows, ncols = Rasters._balance_grid(nplots)
        nrows, ncols, plotinds
    elseif ncols isa Int && nrows isa Int
        nrows, ncols, eachindex(series)
    else
        @error("The provided combination of `ncols::$(typeof(ncols)) and nrows::$(typeof(nrows)) is unsupported. Please either set both to `Makie.automatic` or provide `Int` values.")
    end

    layout = GridLayout(gp, nrows, ncols)

    axes = map(enumerate(plotinds)) do (i, ip)
        ax, plt = Rasters.rplot(layout[fldmod1(i, ncols)...], lift_layer(series, ip); kw...)
        if fldmod1(i, ncols)[1] != nrows
            hidexdecorations!(ax; HIDE_DEC...)
        end
        if fldmod1(i, ncols)[2] != 1
            hideydecorations!(ax; HIDE_DEC...)
        end
        ax
    end
    linkaxes!(axes...)

    return layout
end
function Rasters.rplot(gp::GridPosition, stack::Union{RasterStack,Observable{<:RasterStack}};
    ncols=Makie.automatic,
    nrows=Makie.automatic,
    colormap=nothing,
    colorrange=Makie.automatic,
    link_colorrange=false,
    link_axes=true,
    axis=(;),
    kw...
)
    val_stack = Makie.to_value(stack)
    @assert (length(size(val_stack)) == 2 || size(val_stack, 3) == 1)

    nrows, ncols = if ncols isa Makie.Automatic && nrows isa Makie.Automatic
        Rasters._balance_grid(length(propertynames(val_stack)))
    elseif ncols isa Int && nrows isa Int
        @assert ncols * nrows ≥ length(propertynames(val_stack))
        nrows, ncols
    else
        @error("The provided combination of `ncols::$(typeof(ncols)) and nrows::$(typeof(nrows)) is unsupported. Please either set both to `Makie.automatic` or provide integer values.")
    end

    layout = GridLayout(gp, nrows, ncols)

    axs = [] # avoid ambiguity with Base.axes
    plots = []

    numbered_names = collect(enumerate(propertynames(val_stack)))
    Makie.broadcast_foreach(numbered_names, colormap, colorrange) do (i, name), cmap, crange
        raster = lift_layer(stack, name)
        if length(size(to_value(raster))) == 2
            raster = raster
        elseif length(size(to_value(raster))) == 3
            if size(to_value(raster), 3) == 1
                raster = lift_layer(raster, Rasters.Band(1))
            else
                @error "You can't plot a RasterStack of 3-D rasters using `rplot`. Please provide a stack of 2D rasters instead, or 3D rasters with a singleton third dimension."
            end
        else
            @error "`rplot` cannot plot a Raster of dimension $(size(raster)). Please provide a stack of 2D rasters instead."
        end

        ax, plt = Rasters.rplot(layout[fldmod1(i, ncols)...], raster;
            colormap=cmap, colorrange=crange, axis=axis, kw...
        )

        if fldmod1(i, ncols)[1] != nrows
            hidexdecorations!(ax; HIDE_DEC...)
        end
        if fldmod1(i, ncols)[2] != 1
            hideydecorations!(ax; HIDE_DEC...)
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
function Rasters.rplot(raster::Union{AbstractRaster{T,2},Observable{<:AbstractRaster{T,2}}};
    figure=(;), kw...
) where T
    figure = Figure(; figure...)
    axis, plot = Rasters.rplot(figure[1, 1], raster; kw...)
    return Makie.FigureAxisPlot(figure, axis, plot)
end
function Rasters.rplot(raster::Union{AbstractRaster{T,3},Observable{<:AbstractRaster{T,3}}};
    figure=(;), colormap=nothing, colorrange=Makie.automatic, kw...
) where T
    figure = Figure(; figure...)
    layout = Rasters.rplot(figure[1, 1], raster; colormap, colorrange, kw...)
    # if draw_title
    #     Label(layout[0, 1:Makie.ncols(layout)], raster_title; fontsize = get(figure.scene.attributes, (:Axis :titlesize), 16), font = get(figure.scene.attributes, (:Axis, :titlefont), :bold))
    # end
    return figure
end
function Rasters.rplot(series::Union{AbstractRasterSeries{T},Observable{<:AbstractRasterSeries{T}}};
    figure=(;), colormap=nothing, colorrange=Makie.automatic, kw...
) where T
    figure = Figure(; figure...)
    layout = Rasters.rplot(figure[1, 1], series; colormap, colorrange, kw...)
    return figure
end
function Rasters.rplot(raster::Union{RasterStack,Observable{<:RasterStack}};
    figure=(;), colormap=nothing, colorrange=Makie.automatic, kw...
)
    figure = Figure(; figure...)
    layout = Rasters.rplot(figure[1, 1], raster; colormap, colorrange, kw...)
    # if draw_title
    #     Label(layout[0, 1:Makie.ncols(layout)], raster_title; fontsize = get(figure.scene.attributes, (:Axis :titlesize), 16), font = get(figure.scene.attributes, (:Axis, :titlefont), :bold))
    # end
    return figure
end


# Some small diversions from the DimensionalData.jl recipes
# So 3d plots and series just plot as a single heatmap by default.
# TODO: `plot` should call `rplot`

# 3d rasters are a little more complicated - if dim3 is a singleton, then heatmap, otherwise volume
function Makie.plottype(raster::AbstractRaster{<:Union{Missing,Real},3})
    if size(raster, 3) == 1
        Makie.Heatmap
    else
        Makie.Volume
    end
end

function Makie.convert_arguments(t::Makie.PointBased, A::AbstractRaster{<:Any,1})
    return Makie.convert_arguments(t, _prepare_dimarray(A))
end
function Makie.convert_arguments(t::Makie.PointBased, A::AbstractRaster{<:Number,2})
    return Makie.convert_arguments(t, _prepare_dimarray(A))
end
@static if isdefined(Makie, :SurfaceLike)

    function Makie.convert_arguments(t::SurfaceLike, A::AbstractRaster{var"#s115", 2, D}) where {var"#s115", D<:Tuple}
        return Makie.convert_arguments(t, _prepare_dimarray(A))
    end
else # surfacelike is not a thing
    Makie.convert_arguments(t::Makie.VertexGrid, A::AbstractRaster{<: Any, 2}) = Makie.convert_arguments(t, _prepare_dimarray(A))
    Makie.convert_arguments(t::Makie.CellGrid, A::AbstractRaster{<: Any, 2}) = Makie.convert_arguments(t, _prepare_dimarray(A))
    Makie.convert_arguments(t::Makie.ImageLike, A::AbstractRaster{<: Any, 2}) = Makie.convert_arguments(t, _prepare_dimarray(A))
end

function Makie.convert_arguments(t::Makie.VolumeLike, A::AbstractRaster{<:Any,3}) 
    return Makie.convert_arguments(t, _prepare_dimarray(A))
end
# allow plotting 3d rasters with singleton third dimension (basically 2d rasters)
function Makie.convert_arguments(x::Makie.ConversionTrait, raster::AbstractRaster{<:Union{Real,Missing},3})
    D = _series_dim(raster)
    nplots = size(raster, D)
    if nplots > 1
        # Plot as a RasterSeries
        return Makie.convert_arguments(x, slice(raster, D))
    else
        return Makie.convert_arguments(x, view(raster, rebuild(D, 1)))
    end
end
function Makie.convert_arguments(x::Makie.ConversionTrait, series::AbstractRasterSeries)
    return Makie.convert_arguments(x, first(series))
end

function _series_dim(A)
    spatialdims = (X(), Y(), Z())
    last((dims(A, spatialdims)..., otherdims(A, spatialdims)...))
end

_prepare_dimarray(A) = DimArray(map(x -> _convert_with_missing(x, missingval(A)), A))

_convert_with_missing(x::Real, missingval) = isequal(x, missingval) || ismissing(x) ? NaN32 : Float32(x)
_convert_with_missing(x, missingval) = isequal(x, missingval) ? missing : x
