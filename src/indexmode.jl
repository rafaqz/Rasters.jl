"""
A [`SampledGrid`](@ref) where all cells are the same size and evenly spaced.
These grids will often be paired with a range, but may also be paired with a vector.

## Fields
- `order::Order`: `Order` trait indicating array and index order
- `locus::Locus`: `Locus` trait indicating the position of the indexed point within the cell step
- `span::Number`: the size of a grid step, such as 1u"km" or `Month(1)`
"""
struct ProjectedGrid{O<:Order,Sp,Sa,C,IC} <: AbstractSampledGrid{O,Sp,Sa}
    order::O
    span::Sp
    sampling::Sa
    crs::C
    selectorcrs::IC
end
ProjectedGrid(; order=Ordered(), span=UnknownSpan(), 
              sampling=PointSampling(), crs, selectorcrs=nothing) =
    ProjectedGrid(order, span, sampling, crs, selectorcrs)

crs(grid::ProjectedGrid, dim) = crs(grid)
crs(grid::ProjectedGrid) = grid.crs

selectorcrs(grid::ProjectedGrid, dim) = selectorcrs(grid)
selectorcrs(grid::ProjectedGrid) = grid.selectorcrs

selectorcrs(grid::ProjectedGrid) = grid.selectorcrs

rebuild(g::ProjectedGrid, order=order(g), span=span(g), 
        sampling=sampling(g), crs=crs(g), selectorcrs=selectorcrs(g)) =
    ProjectedGrid(order, span, sampling, crs, selectorcrs)
