"""
    crs(x::Raster)

Get the projected coordinate reference system of a `Y` or `X` `Dimension`,
or of the `Y`/`X` dims of an `AbstractRaster`.

For [`Mapped`](@ref) lookup this may be `nothing` as there may be no projected
coordinate reference system at all.
See [`setcrs`](@ref) to set it manually.
"""
function GeoInterface.crs(obj::Union{<:AbstractRaster,<:AbstractRasterStack,<:AbstractRasterSeries, <:DimTuple})
    each_dim_crs = map(crs, dims(obj))
    firstcrs = findfirst(!isnothing, each_dim_crs)
    if isnothing(firstcrs)
        return nothing
    else
        for (dim, crs) in zip(dims(obj), each_dim_crs)
            if !isnothing(crs) && crs !== each_dim_crs[firstcrs]
                throw(ArgumentError("""
                All dimensions must have the same crs, but dims $(name(dim)) and $(name(dims(obj, firstcrs)))
                have different CRS:
                $(each_dim_crs[firstcrs])

                and

                $(crs)
                """
                ))
            end
        end
        return each_dim_crs[firstcrs]
    end
end
GeoInterface.crs(dim::Dimension) = crs(lookup(dim))

"""
    mappedcrs(x)

Get the mapped coordinate reference system for the `Y`/`X` dims of an array.

In [`Projected`](@ref) lookup this is used to convert [`Selector`]($DDselectordocs)
values form the mappedcrs defined projection to the underlying projection, and to
show plot axes in the mapped projection.

In `Mapped` lookup this is the coordinate reference system of the index values.
See [`setmappedcrs`](@ref) to set it manually.
"""
function mappedcrs end
function mappedcrs(obj)
    if hasdim(obj, Y)
        mappedcrs(dims(obj, Y))
    elseif hasdim(obj, X)
        mappedcrs(dims(obj, X))
    else
        nothing
    end
end
mappedcrs(dim::Dimension) = mappedcrs(lookup(dim))


"""
    setcrs(x, crs)

Set the crs of a `Raster`, `RasterStack`, `Tuple` of `Dimension`, or a `Dimension`.
The `crs` is expected to be a GeoFormatTypes.jl `CRS` or `Mixed` `GeoFormat` type
"""
setcrs(x::Union{<:AbstractRaster,AbstractRasterStack}, crs) = set(x, setcrs(dims(x), crs)...)
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
setmappedcrs(x::Union{<:AbstractRaster,AbstractRasterStack}, mappedcrs) =
    set(x, setmappedcrs(dims(x), mappedcrs)...)
setmappedcrs(dims::DimTuple, mappedcrs) = map(d -> setmappedcrs(d, mappedcrs), dims)
function setmappedcrs(dim::Dimension, mappedcrs)
    rebuild(dim, setmappedcrs(parent(dim), mappedcrs; dim))
end
setmappedcrs(l::AbstractProjected, mappedcrs; dim) = rebuild(l; mappedcrs, dim=basetypeof(dim)())
function setmappedcrs(l::Sampled, mappedcrs; dim)
    dim isa Union{XDim,YDim} ? Mapped(l; mappedcrs, dim) : l
end
setmappedcrs(A::AbstractArray, mappedcrs; dim=nothing) = A

