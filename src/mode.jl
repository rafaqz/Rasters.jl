"""
    Projected(order::Order, span, sampling, crs, usercrs)
    Projected(; order=Ordered(), span=UnknownSpan(), sampling=Points(), crs, usercrs=nothing)

An `AbstractSampled` mode with projections attached.

Fields and behaviours are identical to `Sampled` with the addition of
`crs` and `usercrs` fields.

If both `crs` and `usercrs` fields contain crs data (in a GeoFormat wrapper 
from GeoFormatTypes.jl) the selector inputs and plot axes will be converted 
from and to the specified `usercrs` projection automatically. A common use case
would be to pass `usercrs=EPSG(4326)` to the constructor when loading eg. a GDALarray:

```julia
GDALarray(filename; usercrs=EPSG(4326))
```

The underlying `crs` will be detected by GDAL.

If `usercrs` is not supplied (ie. `isa Nothing`), the base index will be shown on plots, 
and selectors will need to use whatever format it is in.
"""
struct Projected{O<:Order,Sp,Sa,C,IC} <: AbstractSampled{O,Sp,Sa}
    order::O
    span::Sp
    sampling::Sa
    crs::C
    usercrs::IC
end
Projected(; order=Ordered(), span=UnknownSpan(), 
          sampling=Points(), crs, usercrs=nothing) =
    Projected(order, span, sampling, crs, usercrs)

crs(mode::Projected, dim) = crs(mode)
crs(mode::Projected) = mode.crs

usercrs(mode::Projected, dim) = usercrs(mode)
usercrs(mode::Projected) = mode.usercrs

rebuild(g::Projected, order=order(g), span=span(g), 
        sampling=sampling(g), crs=crs(g), usercrs=usercrs(g)) =
    Projected(order, span, sampling, crs, usercrs)

"""
    Converted(order::Order, span, sampling, crs, dimcrs)
    Converted(; order=Ordered(), span=UnknownSpan(), sampling=Points(), crs, dimcrs)

An `AbstractSampled` mode with projections, where the dimension has already been converted
to another projection as a vector, usually `EPSG(4326)`.

Fields and behaviours are identical to `Sampled` with the addition of
`crs` and `dimcrs` fields.

The dimension will be indexed as for `Sampled`, but to save in another format the
underlying projection will be used.

```julia
GDALarray(filename; usercrs=EPSG(4326))
```

The underlying `crs` will be detected by GDAL.

If `usercrs` is not supplied (ie. `isa Nothing`), the base index will be shown on plots, 
and selectors will need to use whatever format it is in.
"""
struct Converted{O<:Order,Sp,Sa,C,DC} <: AbstractSampled{O,Sp,Sa}
    order::O
    span::Sp
    sampling::Sa
    crs::C
    dimcrs::DC
end
Converted(; order=Ordered(), span=UnknownSpan(), 
          sampling=Points(), crs, dimcrs) =
    Converted(order, span, sampling, crs, dimcrs)

crs(mode::Converted, dim) = crs(mode)
crs(mode::Converted) = mode.crs

dimcrs(mode::Converted, dim) = dimcrs(mode)
dimcrs(mode::Converted) = mode.dimcrs

rebuild(g::Converted, order=order(g), span=span(g), 
        sampling=sampling(g), crs=crs(g), dimcrs=dimcrs(g)) =
    Converted(order, span, sampling, crs, dimcrs)

"""
    LatLon(order, span, sampling)
    LatLon(; order=Ordered(), span=UnknownSpan(), sampling=Points())

An `AbstractSampled` mode for standard latitude/longitude dimensions.
"""
struct LatLon{O<:Order,Sp,Sa} <: AbstractSampled{O,Sp,Sa}
    order::O
    span::Sp
    sampling::Sa
end
LatLon(; order=Ordered(), span=UnknownSpan(), sampling=Points()) =
    LatLon(order, span, sampling)

crs(mode::LatLon, args...) = EPSG(4326)

rebuild(g::LatLon, order=order(g), span=span(g), sampling=sampling(g)) =
    LatLon(order, span, sampling)
