# These methods are only available if ArchGDAL is loaded.
# Otherwise ProjectedGrid selector crs field must be `nothing`,
# and no reprojection can occur.

sel2indices(grid::ProjectedGrid, dim::Dimension, sel::Contains{<:Number}) = begin
    selval = reproject(grid, dim, val(sel))
    DD.contains(dim, rebuild(sel, selval))
end
sel2indices(grid::ProjectedGrid, dim::Dimension, sel::At{<:Number}) = begin
    selval = reproject(grid, dim, val(sel))
    DD.at(dim, rebuild(sel, selval))
end
sel2indices(grid::ProjectedGrid, dim::Dimension, sel::Between) = begin
    selval = map(v -> reproject(grid, dim, v), val(sel))
    DD.between(dim, rebuild(sel, selval))
end

# TODO pass whole vector to ArchGDAL
reproject(grid::ProjectedGrid, dim::Dimension, val::Number) = 
    reproject(grid, selectorcrs(grid), crs(grid), dim, val::Number) 
reproject(grid::ProjectedGrid, target::Nothing, source, dim, selval::Number) = selval 
reproject(grid::ProjectedGrid, target::GeoFormat, source::GeoFormat, dim::Lat, selval::Number) = 
    ArchGDAL.reproject([(selval, zero(selval))], target, source)[1][2]
reproject(grid::ProjectedGrid, target::GeoFormat, source::GeoFormat, dim::Lon, selval::Number) = 
    ArchGDAL.reproject([(zero(selval), selval)], target, source)[1][1]
