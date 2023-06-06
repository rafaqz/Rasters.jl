"""
    AbstractProjected <: AbstractSampled

Abstract supertype for projected index lookups.
"""
abstract type AbstractProjected{T,O,Sp,Sa} <: AbstractSampled{T,O,Sp,Sa} end

struct AutoDim end

# For now we just remove CRS on GPU - it often contains strings
function Adapt.adapt_structure(to, l::AbstractProjected)
    sampled = Sampled(parent(l), order(l), span(l), sampling(l), metadata(l))
    return Adapt.adapt_structure(to, sampled)
end

crs(lookup::LookupArray) = nothing
mappedcrs(lookup::LookupArray) = nothing

# When the lookup is formatted with an array we match the `dim` field with the
# wrapper dimension. We will need this later for e.g. projecting with GDAL,
# where we need to know which lookup is X and which is Y
function Dimensions.format(m::AbstractProjected, D::Type, index, axis::AbstractRange)
    i = Dimensions._format(index, axis)
    o = Dimensions._format(order(m), D, index)
    sp = Dimensions._format(span(m), D, index)
    sa = Dimensions._format(sampling(m), sp, D, index)
    dim = Dimensions.basetypeof(D)()
    x = rebuild(m; data=i, order=o, span=sp, sampling=sa, dim=dim)
    return x
end

"""
    Projected <: AbstractProjected

    Projected(order, span, sampling, crs, mappedcrs)
    Projected(; order=AutoOrder(), span=AutoSpan(), sampling=AutoSampling(), crs, mappedcrs=nothing)

An [`AbstractSampled`]($DDabssampleddocs) `LookupArray` with projections attached.

Fields and behaviours are identical to [`Sampled`]($DDsampleddocs)
with the addition of `crs` and `mappedcrs` fields.

If both `crs` and `mappedcrs` fields contain CRS data (in a `GeoFormat` wrapper
from GeoFormatTypes.jl) the selector inputs and plot axes will be converted
from and to the specified `mappedcrs` projection automatically. A common use case
would be to pass `mappedcrs=EPSG(4326)` to the constructor when loading eg. a GDALarray:

```julia
GDALarray(filename; mappedcrs=EPSG(4326))
```

The underlying `crs` will be detected by GDAL.

If `mappedcrs` is not supplied (ie. `mappedcrs=nothing`), the base index will be
shown on plots, and selectors will need to use whatever format it is in.
"""
struct Projected{T,A<:AbstractVector{T},O<:Order,Sp<:Span,Sa<:Sampling,MD,PC,MC,D} <: AbstractProjected{T,O,Sp,Sa}
    data::A
    order::O
    span::Sp
    sampling::Sa
    metadata::MD
    crs::PC
    mappedcrs::MC
    dim::D
end
function Projected(data=AutoIndex();
    order=AutoOrder(), span=AutoSpan(), sampling=AutoSampling(),
    metadata=NoMetadata(), crs, mappedcrs=nothing, dim=AutoDim()
)
    Projected(data, order, span, sampling, metadata, crs, mappedcrs, dim)
end
function Projected(l::Sampled;
    order=order(l), span=span(l), sampling=sampling(l),
    metadata=metadata(l), crs, mappedcrs=nothing, dim=AutoDim()
)
    Projected(parent(l), order, span, sampling, metadata, crs, mappedcrs, dim)
end

crs(lookup::Projected) = lookup.crs
mappedcrs(lookup::Projected) = lookup.mappedcrs
dim(lookup::Projected) = lookup.dim

function LA.selectindices(l::Projected, sel::Contains)
    selval = reproject(mappedcrs(l), crs(l), dim(l), val(sel))
    LA.contains(l, rebuild(sel; val=selval))
end
function LA.selectindices(l::Projected, sel::At)
    selval = reproject(mappedcrs(l), crs(l), dim(l), val(sel))
    LA.at(l, rebuild(sel; val=selval))
end
function LA.selectindices(l::Projected, sel::Between)
    selval = map(v -> reproject(mappedcrs(l), crs(l), dim(l), v), val(sel))
    LA.between(l, rebuild(sel; val=selval))
end

"""
    Mapped <: AbstractProjected

    Mapped(order, span, sampling, crs, mappedcrs)
    Mapped(; order=AutoOrder(), span=AutoSpan(), sampling=AutoSampling(), crs=nothing, mappedcrs)

An [`AbstractSampled`]($DDabssampleddocs) `LookupArray`, where the dimension index has
been mapped to another projection, usually lat/lon or `EPSG(4326)`.
`Mapped` matches the dimension format commonly used in netcdf files.

Fields and behaviours are identical to [`Sampled`]($DDsampleddocs) with the addition of
`crs` and `mappedcrs` fields.

The mapped dimension index will be used as for [`Sampled`]($DDsampleddocs),
but to save in another format the underlying `crs` may be used to convert it.
"""
struct Mapped{T,A<:AbstractVector{T},O<:Order,Sp<:Span,Sa<:Sampling,MD,PC,MC,D} <: AbstractProjected{T,O,Sp,Sa}
    data::A
    order::O
    span::Sp
    sampling::Sa
    metadata::MD
    crs::PC
    mappedcrs::MC
    dim::D
end
function Mapped(data=AutoIndex();
    order=AutoOrder(), span=AutoSpan(), sampling=AutoSampling(),
    metadata=NoMetadata(), crs=nothing, mappedcrs, dim=AutoDim()
)
    Mapped(data, order, span, sampling, metadata, crs, mappedcrs, dim)
end
function Mapped(l::Sampled;
    order=order(l), span=span(l), sampling=sampling(l),
    metadata=metadata(l), crs=nothing, mappedcrs, dim=AutoDim()
)
    Mapped(parent(l), order, span, sampling, metadata, crs, mappedcrs, dim)
end

crs(lookup::Mapped) = lookup.crs
mappedcrs(lookup::Mapped) = lookup.mappedcrs
dim(lookup::Mapped) = lookup.dim

"""
    convertlookup(dstlookup::Type{<:LookupArray}, x)

Convert the dimension lookup between `Projected` and `Mapped`.
Other dimension lookups pass through unchanged.

This is used to e.g. save a netcdf file to GeoTiff.
"""
convertlookup(T::Type{<:LookupArray}, A::AbstractDimArray) =
    rebuild(A; dims=convertlookup(T, dims(A)))
convertlookup(T::Type{<:LookupArray}, dims::Tuple) = map(d -> convertlookup(T, d), dims)
convertlookup(T::Type{<:LookupArray}, d::Dimension) = rebuild(d, convertlookup(T, lookup(d)))
# Non-projected LookupArray lookupss pass through
convertlookup(::Type, lookup::LookupArray) = lookup
# AbstractProjected passes through if it's the same as dstlookup
convertlookup(::Type{T1}, lookup::T2) where {T1,T2<:T1} = lookup
# Otherwise AbstractProjected needs ArchGDAL
function convertlookup(::Type{<:Mapped}, l::Projected)
    newindex = reproject(crs(l), mappedcrs(l), dim(l), index(l))
    # We use Explicit mode and make a bounds matrix
    # This way the bounds can be saved correctly to NetCDF
    newbounds = reproject(crs(l), mappedcrs(l), dim(l), Dimensions.dim2boundsmatrix(l))
    return Mapped(newindex,
        order=order(l),
        span=Explicit(newbounds),
        sampling=sampling(l),
        metadata=metadata(l),
        crs=crs(l),
        mappedcrs=mappedcrs(l),
        dim=dim(l),
    )
end
function convertlookup(::Type{<:Projected}, l::Mapped)
    newindex = _projectedrange(l)
    return Projected(newindex;
        order=order(l),
        span=Regular(step(newindex)),
        sampling=sampling(l),
        metadata=metadata(l),
        crs=crs(l),
        mappedcrs=mappedcrs(l),
        dim=dim(l),
    )
end



_projectedrange(l::Projected) = LinRange(first(l), last(l), length(l))
_projectedrange(l::Mapped) = _projectedrange(span(l), crs(l), l)
function _projectedrange(span, crs, l::Mapped)
    start, stop = reproject(mappedcrs(l), crs, dim(l), [first(l), last(l)])
    LinRange(start, stop, length(l))
end
_projectedrange(::Regular, crs::Nothing, l::Mapped) = LinRange(first(l), last(l), length(l))
function _projectedrange(::T, crs::Nothing, l::Mapped) where T<:Union{Irregular,Explicit}
    error("Cannot convert a Mapped $T index to Projected when crs is nothing")
end

"""
    setcrs(x, crs)

Set the crs of a `Raster`, `RasterStack`, `Tuple` of `Dimension`, or a `Dimension`.
The `crs` is expected to be a GeoFormatTypes.jl `CRS` or `Mixed` `GeoFormat` type
"""
setcrs(dims::DimTuple, crs) = map(d -> setcrs(d, crs), dims)
function setcrs(dim::Dimension, crs)
    rebuild(dim, setcrs(parent(dim), crs; dim=basetypeof(dim)()))
end
setcrs(l::AbstractProjected, crs; dim=nothing) = rebuild(l; crs)
function setcrs(l::Sampled, crs; dim)
    dim isa Union{XDim,YDim} ? Projected(l; crs, dim) : l
end
setcrs(A::AbstractArray, crs; dim=nothing) = A

"""
    setmappedcrs(x, crs)

Set the mapped crs of a `Raster`, a `RasterStack`, a `Tuple`
of `Dimension`, or a `Dimension`.
The `crs` is expected to be a GeoFormatTypes.jl `CRS` or `Mixed` `GeoFormat` type
"""
setmappedcrs(dims::DimTuple, mappedcrs) = map(d -> setmappedcrs(d, mappedcrs), dims)
function setmappedcrs(dim::Dimension, mappedcrs)
    rebuild(dim, setmappedcrs(parent(dim), mappedcrs; dim))
end
setmappedcrs(l::AbstractProjected, mappedcrs; dim) = rebuild(l; mappedcrs, dim=basetypeof(dim)())
function setmappedcrs(l::Sampled, mappedcrs; dim)
    dim isa Union{XDim,YDim} ? Mapped(l; mappedcrs, dim) : l
end
setmappedcrs(A::AbstractArray, mappedcrs; dim=nothing) = A


"""
    mappedbounds(x)

Get the bounds converted to the [`mappedcrs`](@ref) value.

Whithout ArchGDAL loaded, this is just the regular bounds.
"""
function mappedbounds end

mappedbounds(dims::Tuple) = map(mappedbounds, dims)
mappedbounds(dim::Dimension) = mappedbounds(parent(dim), dim)
mappedbounds(::LookupArray, dim) = bounds(dim)
mappedbounds(lookup::Projected, dim) = mappedbounds(mappedcrs(lookup), lookup, dim)
mappedbounds(mappedcrs::Nothing, lookup::Projected, dim) =
    error("No mappedcrs attached to $(name(dim)) dimension")
mappedbounds(mappedcrs::GeoFormat, lookup::Projected, dim) =
    _sort(reproject(crs(lookup), mappedcrs, dim, bounds(dim)))

projectedbounds(dims::Tuple) = map(projectedbounds, dims)
projectedbounds(dim::Dimension) = projectedbounds(parent(dim), dim)
projectedbounds(::LookupArray, dim) = bounds(dim)
projectedbounds(lookup::Mapped, dim) = projectedbounds(crs(lookup), lookup, dim)
projectedbounds(crs::Nothing, lookup::Mapped, dim) =
    error("No projection crs attached to $(name(dim)) dimension")
projectedbounds(crs::GeoFormat, lookup::Mapped, dim) =
    _sort(reproject(mappedcrs(lookup), crs, dim, bounds(dim)))

_sort((a, b)) = a <= b ? (a, b) : (b, a)

"""
    mappedindex(x)

Get the index value of a dimension converted to the `mappedcrs` value.

Whithout ArchGDAL loaded, this is just the regular dim value.
"""
function mappedindex end

mappedindex(dims::Tuple) = map(mappedindex, dims)
mappedindex(dim::Dimension) = _mappedindex(parent(dim), dim)

_mappedindex(::LookupArray, dim::Dimension) = index(dim)
_mappedindex(lookup::Projected, dim::Dimension) = _mappedindex(mappedcrs(lookup), lookup, dim)
_mappedindex(mappedcrs::Nothing, lookup::Projected, dim) =
    error("No mappedcrs attached to $(name(dim)) dimension")
_mappedindex(mappedcrs::GeoFormat, lookup::Projected, dim) =
    reproject(crs(dim), mappedcrs, dim, index(dim))

projectedindex(dims::Tuple) = map(projectedindex, dims)
projectedindex(dim::Dimension) = _projectedindex(parent(dim), dim)

_projectedindex(::LookupArray, dim::Dimension) = index(dim)
_projectedindex(lookup::Mapped, dim::Dimension) = _projectedindex(crs(lookup), lookup, dim)
_projectedindex(crs::Nothing, lookup::Mapped, dim::Dimension) =
    error("No projection crs attached to $(name(dim)) dimension")
_projectedindex(crs::GeoFormat, lookup::Mapped, dim::Dimension) =
    reproject(mappedcrs(dim), crs, dim, index(dim))
