# These methods are only available if ArchGDAL is loaded.
# Otherwise ProjectedIndex selector crs field must be `nothing`,
# and no reprojection can occur.

sel2indices(mode::ProjectedIndex, dim::Dimension, sel::Contains{<:Number}) = begin
    selval = reproject(mode, dim, val(sel))
    println((name(dim), selval))
    DD.contains(dim, rebuild(sel, selval))
end
sel2indices(mode::ProjectedIndex, dim::Dimension, sel::At{<:Number}) = begin
    selval = reproject(mode, dim, val(sel))
    println((name(dim), selval))
    DD.at(dim, rebuild(sel, selval))
end
sel2indices(mode::ProjectedIndex, dim::Dimension, sel::Between) = begin
    selval = map(v -> reproject(mode, dim, v), val(sel))
    println((name(dim), selval))
    DD.between(dim, rebuild(sel, selval))
end

# TODO pass whole vector to ArchGDAL
reproject(mode::ProjectedIndex, dim::Dimension, val::Number) =
    reproject(mode, usercrs(mode), crs(mode), dim, val::Number)
reproject(mode::ProjectedIndex, target::Nothing, source, dim, selval::Number) = selval
reproject(mode::ProjectedIndex, target::GeoFormat, source::GeoFormat, dim::Lat, selval::Number) =
    ArchGDAL.reproject([(selval, zero(selval))], target, source)[1][2]
reproject(mode::ProjectedIndex, target::GeoFormat, source::GeoFormat, dim::Lon, selval::Number) =
    ArchGDAL.reproject([(zero(selval), selval)], target, source)[1][1]
