export reproject, convertmode, userval, userbounds

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
    newval = reproject(crs(m), usercrs(m), dim, val(dim))
    newbounds = reproject(crs(m), usercrs(m), dim, bounds(dim))
    newmode = Mapped(
        order=order(m),
        span=Irregular(newbounds),
        sampling=sampling(m),
        projectedcrs=crs(m),
        mappedcrs=usercrs(m),
    )
    rebuild(dim; val=newval, mode=newmode)
end
convertmode(dstmode::Type{Projected}, srcmode::Type{Mapped}, dim::Dimension) where M = begin
    m = mode(dim)
    start, stop = reproject(mappedcrs(m), projectedcrs(m), dim, [first(dim), last(dim)]) newval = LinRange(start, stop, length(dim))
    newmode = Projected(
        order=order(m),
        span=Regular(step(newval)),
        sampling=sampling(m),
        crs=projectedcrs(m),
        usercrs=mappedcrs(m),
    )
    rebuild(dim; val=newval, mode=newmode)
end


# Add Lat/Lon methods to reproject bounds and val
userbounds(dim::Union{Lat,Lon}) = userbounds(mode(dim), dim)
userbounds(::Projected, dim::Union{Lat,Lon}) =
    reproject(crs(dim), usercrs(dim), dim, bounds(dim))
userbounds(::IndexMode, dim::Union{Lat,Lon}) = bounds(dim)

userval(dim::Union{Lat,Lon}) = userval(mode(dim), dim)
userval(::Projected, dim::Union{Lat,Lon}) =
    reproject(crs(dim), usercrs(dim), dim, val(dim))
userval(::IndexMode, dim::Union{Lat,Lon}) = val(dim)
