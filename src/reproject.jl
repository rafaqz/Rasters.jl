export reproject, convertmode


# These methods are only available if ArchGDAL is loaded.
# Otherwise Projected selector crs field must be `nothing`,
# and no reprojection can occur.

sel2indices(mode::Projected, dim::Dimension, sel::Contains{<:Number}) = begin
    selval = reproject(usercrs(mode), crs(mode), dim, val(sel))
    DD.contains(dim, rebuild(sel, selval))
end
sel2indices(mode::Projected, dim::Dimension, sel::At{<:Number}) = begin
    selval = reproject(usercrs(mode), crs(mode), dim, val(sel))
    DD.at(dim, rebuild(sel, selval))
end
sel2indices(mode::Projected, dim::Dimension, sel::Between) = begin
    selval = map(v -> reproject(usercrs(mode), crs(mode), dim, v), val(sel))
    DD.between(dim, rebuild(sel, selval))
end

"""
`reproject` uses ArchGDAL.reproject, but implemented for a reprojecting
a single dimension at a time.
"""
reproject(source, target, dim::Dimension, val) = val
reproject(source::GeoFormat, target::GeoFormat, dim::Lon, val::Number) =
    ArchGDAL.reproject((Float64(val), 0.0), source, target; order=:trad)[1]
reproject(source::GeoFormat, target::GeoFormat, dim::Lat, val::Number) =
    ArchGDAL.reproject((0.0, Float64(val)), source, target; order=:trad)[2]

reproject(source::GeoFormat, target::GeoFormat, ::Lon, vals::AbstractArray) =
    [r[1] for r in ArchGDAL.reproject([(Float64(v), 0.0) for v in vals], source, target; order=:trad)]
reproject(source::GeoFormat, target::GeoFormat, dim::Lat, vals::AbstractArray) =
    [r[2] for r in ArchGDAL.reproject([(0.0, Float64(v)) for v in vals], source, target; order=:trad)]

reproject(source::GeoFormat, target::GeoFormat, ::Lon, vals::Tuple) =
    Tuple(r[1] for r in ArchGDAL.reproject([(Float64(v), 0.0) for v in vals], source, target; order=:trad))
reproject(source::GeoFormat, target::GeoFormat, dim::Lat, vals::Tuple) =
    Tuple(r[2] for r in ArchGDAL.reproject([(0.0, Float64(v)) for v in vals], source, target; order=:trad))


convertmode(dstmode::Type{Mapped}, srcmode::Type{Projected}, dim::Dimension) where M = begin
    m = mode(dim)
    newindex = reproject(projectedcrs(m), mappedcrs(m), dim, val(dim))
    newbounds = reproject(projectedcrs(m), mappedcrs(m), dim, bounds(dim))
    newmode = Mapped(
        order=order(m),
        span=Irregular(newbounds),
        sampling=sampling(m),
        projectedcrs=projectedcrs(m),
        mappedcrs=mappedcrs(m),
    )
    rebuild(dim; val=newindex, mode=newmode)
end
convertmode(dstmode::Type{Projected}, srcmode::Type{Mapped}, dim::Dimension) where M = begin
    m = mode(dim)
    newindex = _projectedrange(m, dim)
    newmode = Projected(
        order=order(m),
        span=Regular(step(newindex)), sampling=sampling(m),
        projectedcrs=projectedcrs(m),
        mappedcrs=mappedcrs(m),
    )
    rebuild(dim; val=newval, mode=newmode)
end

_projectedrange(::Projected, dim) = LinRange(first(dim), last(dim), length(dim))
_projectedrange(m::Mapped, dim) = 
    _projectedrange(span(m), projectedcrs(m), m, dim) 
_projectedrange(span, projectedcrs, m::Mapped, dim) = begin
    start, stop = reproject(mappedcrs(m), projectedcrs, dim, [first(dim), last(dim)]) 
    LinRange(start, stop, length(dim))
end
_projectedrange(::Regular, projectedcrs::Nothing, ::Mapped, dim) =
    LinRange(first(dim), last(dim), length(dim))
_projectedrange(::Irregular, projectedcrs::Nothing, ::Mapped, dim) = 
    error("Cannot convert an Mapped Irregular index to Projected when projectioncrs is nothing")

# Add Lat/Lon methods to reproject bounds and val
projectedbounds(dim::Union{Lat,Lon}) = projectedbounds(mode(dim), dim)
projectedbounds(::Projected, dim::Union{Lat,Lon}) =
    reproject(crs(dim), projectedcrs(dim), dim, bounds(dim))
projectedbounds(::IndexMode, dim::Union{Lat,Lon}) = bounds(dim)

projectedval(dim::Union{Lat,Lon}) = projectedval(mode(dim), dim)
projectedval(::Projected, dim::Union{Lat,Lon}) =
    reproject(crs(dim), projectedcrs(dim), dim, val(dim))
projectedval(::IndexMode, dim::Union{Lat,Lon}) = val(dim)


# Add Lat/Lon methods to reproject bounds and val
mappedbounds(dim::Union{Lat,Lon}) = mappedbounds(mode(dim), dim)
mappedbounds(::Projected, dim::Union{Lat,Lon}) =
    reproject(crs(dim), mappedcrs(dim), dim, bounds(dim))
mappedbounds(::IndexMode, dim::Union{Lat,Lon}) = bounds(dim)

mappedval(dim::Union{Lat,Lon}) = mappedval(mode(dim), dim)
mappedval(::Projected, dim::Union{Lat,Lon}) =
    reproject(crs(dim), mappedcrs(dim), dim, val(dim))
mappedval(::IndexMode, dim::Union{Lat,Lon}) = val(dim)
