abstract type AbstractProjected{O,Sp,Sa} <: AbstractSampled{O,Sp,Sa} end

"""
    Projected(order::Order, span, sampling, projectedcrs, mappedcrs)
    Projected(; order=Ordered(), span=AutoSpan(), sampling=Points(), projectedcrs, mappedcrs=nothing)

An [`AbstractSampled`]($DDabssampleddocs) `IndexMode` with projections attached.

Fields and behaviours are identical to [`Sampled`]($DDsampleddocs)
with the addition of `projectedcrs` and `mappedcrs` fields.

If both `projectedcrs` and `mappedcrs` fields contain CRS data (in a `GeoFormat` wrapper
from GeoFormatTypes.jl) the selector inputs and plot axes will be converted
from and to the specified `mappedcrs` projection automatically. A common use case
would be to pass `mappedcrs=EPSG(4326)` to the constructor when loading eg. a GDALarray:

```julia
GDALarray(filename; mappedcrs=EPSG(4326))
```

The underlying `projectedcrs` will be detected by GDAL.

If `mappedcrs` is not supplied (ie. `mappedcrs=nothing`), the base index will be
shown on plots, and selectors will need to use whatever format it is in.
"""
struct Projected{O<:Order,Sp<:Regular,Sa<:Sampling,PC,MC} <: AbstractProjected{O,Sp,Sa}
    order::O
    span::Sp
    sampling::Sa
    projectedcrs::PC
    mappedcrs::MC
end
Projected(; order=Ordered(), span=Regular(),
          sampling=Points(), projectedcrs, mappedcrs=nothing) =
    Projected(order, span, sampling, projectedcrs, mappedcrs)

projectedcrs(mode::Projected, dim) = projectedcrs(mode)
projectedcrs(mode::Projected) = mode.projectedcrs

mappedcrs(mode::Projected, dim) = mappedcrs(mode)
mappedcrs(mode::Projected) = mode.mappedcrs

rebuild(g::Projected, order=order(g), span=span(g),
        sampling=sampling(g), projectedcrs=projectedcrs(g), mappedcrs=mappedcrs(g)) =
    Projected(order, span, sampling, projectedcrs, mappedcrs)

"""
    Mapped(order::Order, span, sampling, projectedcrs, mappedcrs)
    Mapped(; order=Ordered(), span=AutoSpan(), sampling=Points(), projectedcrs=nothing, mappedcrs)

An [`AbstractSampled`]($DDabssampleddocs) `IndexMode`, where the dimension index has  
been mapped to another projection, usually lat/lon or `EPSG(4326)`.

Fields and behaviours are identical to [`Sampled`]($DDsampleddocs) with the addition of
`projectedcrs` and `mappedcrs` fields.

The mapped dimension index will be used as for [`Sampled`]($DDsampleddocs), 
but to save in another format the underlying `projectioncrs` may be used.
"""
struct Mapped{O<:Order,Sp<:Span,Sa<:Sampling,PC,MC} <: AbstractProjected{O,Sp,Sa}
    order::O
    span::Sp
    sampling::Sa
    projectedcrs::PC
    mappedcrs::MC
end
Mapped(; order=Ordered(), span=AutoSpan(), sampling=Points(), 
       projectedcrs=nothing, mappedcrs) =
    Mapped(order, span, sampling, projectedcrs, mappedcrs)

projectedcrs(mode::Mapped, dim) = projectedcrs(mode)
projectedcrs(mode::Mapped) = mode.projectedcrs

mappedcrs(mode::Mapped, dim) = mappedcrs(mode)
mappedcrs(mode::Mapped) = mode.mappedcrs

rebuild(g::Mapped, order=order(g), span=span(g),
        sampling=sampling(g), projectedcrs=projectedcrs(g), mappedcrs=mappedcrs(g)) =
    Mapped(order, span, sampling, projectedcrs, mappedcrs)

"""
    convertmode(dstmode::Type{<:IndexMode}, x)

Convert the dimension mode between `Projected` and `Mapped`.
Other dimension modes pass through unchanged.

This is used to e.g. save a netcdf file to GeoTiff.
"""
convertmode(dstmode::Type{<:IndexMode}, A::AbstractArray) =
    rebuild(A, data(A), convertmode(dstmode, dims(A)))
convertmode(dstmode::Type{<:IndexMode}, dims::Tuple) =
    map(d -> convertmode(dstmode, d), dims)
convertmode(dstmode::Type{<:IndexMode}, dim::Dimension) =
    convertmode(dstmode, basetypeof(mode(dim)), dim)
# Regular modes pass through
convertmode(dstmode::Type, srcmode::Type{<:IndexMode}, dim::Dimension) = dim
# AbstractProjected passes through if it's the same as dstmode
convertmode(dstmode::Type{M}, srcmode::Type{M}, dim::Dimension) where M<:AbstractProjected = dim
# Otherwise AbstractProjected needs ArchGDAL
convertmode(dstmode::Type, srcmode::Type{<:AbstractProjected}, dim::Dimension) =
    error("Load ArchGDAL.jl to convert projected dimensions")
