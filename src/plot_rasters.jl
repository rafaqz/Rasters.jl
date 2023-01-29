function _balance_grid(nplots)
    ncols = (nplots - 1) ÷ ceil(Int, sqrt(nplots)) + 1
    nrows = (nplots - 1) ÷ ncols + 1
    return ncols, nrows
end

theme_rasters = Theme(
    resolution=(650, 400),
    fontsize=16,
    font="CMU Serif",
    colormap=:plasma, # plasma, imola, batlow
    Axis=(;
        aspect=DataAspect(),
        xlabelsize=18,
        ylabelsize=18,
        xgridstyle=:dash,
        ygridstyle=:dash,
        xtickalign=1,
        ytickalign=1,
        yticksize=10,
        xticksize=10,
        xlabelpadding=-5,
        xlabel="x",
        ylabel="y"),
    Legend=(;
        framecolor=(:orange, 0.5),
        bgcolor=(:white, 0.5)
        ),
    Colorbar=(;
        label ="z",
        ticksize=16,
        tickalign=1,
        spinewidth=0.5),
)
set_theme!(theme_rasters)

function rplot(ras::AbstractRaster{T,2,<:Tuple{D1,D2}};
        X=X, Y=Y,
        xlabel = "lon",
        ylabel = "lat",
        clabel = "z",
        colormap = :plasma,
        nan_color = (:brown, 0.02),
        replace=false
        ) where {T,D1<:SpatialDim,D2<:SpatialDim}
    if replace
        ras = replace_missing(ras, missingval=NaN)
    end

    x, y = lookup(ras, X), lookup(ras, Y)
    d = ras.data
    s = size(d)
    aspect = s[1]/s[2]
    fig = Figure()
    ax = Axis(fig[1,1]; xlabel, ylabel)
    obj = heatmap!(ax, x, y, d; colormap, nan_color)
    Colorbar(fig[1,2], obj; label = clabel)
    colsize!(fig.layout, 1, Aspect(1, aspect))
    colgap!(fig.layout, 10)
    return fig
end

function rplot(ras::Raster{T,3,<:Tuple{<:SpatialDim,<:SpatialDim,D}};
    X=X, Y=Y,
    xlabel = "lon",
    ylabel = "lat",
    clabel = "z",
    colormap = :plasma,
    nan_color = (:brown, 0.02),
    base_sizex = 450,
    base_sizey = 400,
    replace=true
    ) where {T,D}
    if replace
        ras = replace_missing(ras, missingval=NaN)
    end

    x, y = lookup(ras, X), lookup(ras, Y)
    d = ras.data
    s = size(d)
    if s[3]==1
        aspect = s[1]/s[2]

        fig = Figure()
        ax = Axis(fig[1,1]; xlabel, ylabel)
        obj = heatmap!(ax, x, y, d[:,:,1]; colormap, nan_color)
        Colorbar(fig[1,2], obj; label = clabel)
        colsize!(fig.layout, 1, Aspect(1, aspect))
        colgap!(fig.layout, 10)
        return fig
    elseif s[3]>1
        nfacets = s[3]
        ncols, nrows = Rasters._balance_grid(nfacets)
    
        x = lookup(ras, X) 
        y = lookup(ras, Y)
        aspect = s[1]/s[2]

        fig = Figure(resolution=(ncols*base_sizex,nrows*base_sizey))
        gs = [fig[i,j] = GridLayout() for i in 1:nrows for j in 1:ncols]
        for k in 1:nfacets
            ax = Axis(gs[k][1,1]; aspect = 1, title = "$(k)")
            obj = heatmap!(ax, x, y, ras.data[:,:,k]; colormap, nan_color)
            Colorbar(gs[k][1,2], obj)
            colsize!(gs[k], 1, Aspect(1, aspect))
            colgap!(gs[k], 5)
        end
        Label(fig[0,:], string(ras.name), fontsize = 24)
        rowgap!(fig.layout, 10)
        fig
    end
end

export rplot, theme_rasters