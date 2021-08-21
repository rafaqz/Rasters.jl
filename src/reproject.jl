function DD._sel2indices(mode::Projected, dim::Dimension, sel::Contains)
    selval = reproject(mappedcrs(mode), crs(mode), dim, val(sel))
    DD.contains(dim, rebuild(sel, selval))
end
function DD._sel2indices(mode::Projected, dim::Dimension, sel::At)
    selval = reproject(mappedcrs(mode), crs(mode), dim, val(sel))
    DD.at(dim, rebuild(sel, selval))
end
function DD._sel2indices(mode::Projected, dim::Dimension, sel::Between)
    selval = map(v -> reproject(mappedcrs(mode), crs(mode), dim, v), val(sel))
    DD.between(dim, rebuild(sel, selval))
end

"""
    reproject(source::GeoFormat, target::GeoFormat, dim::Dimension, val)

`reproject` uses ArchGDAL.reproject, but implemented for a reprojecting
a single dimension at a time.
"""
reproject(source, target, dim::Dimension, val) = val
function reproject(source::GeoFormat, target::GeoFormat, dim::X, val::Number)
    AG.reproject((val, zero(val)), source, target; order=:trad)[1]
end
function reproject(source::GeoFormat, target::GeoFormat, dim::Y, val::Number)
    AG.reproject((zero(val), val), source, target; order=:trad)[2]
end
function reproject(source::GeoFormat, target::GeoFormat, ::X, vals::AbstractArray) rep = AG.reproject(map(v -> (v, zero(v)), vals), source, target; order=:trad)
    map(r -> r[1], rep)
end
function reproject(source::GeoFormat, target::GeoFormat, dim::Y, vals::AbstractArray)
    rep = AG.reproject(map(v -> (zero(v), v), vals), source, target; order=:trad)
    map(r -> r[2], rep)
end
function reproject(source::GeoFormat, target::GeoFormat, ::X, vals::Tuple)
    reps = AG.reproject([(v, zero(v)) for v in vals], source, target; order=:trad)
    Tuple(r[1] for r in reps)
end
function reproject(source::GeoFormat, target::GeoFormat, dim::Y, vals::Tuple)
    reps = AG.reproject([(zero(v), v) for v in vals], source, target; order=:trad)
    Tuple(r[2] for r in reps)
end
