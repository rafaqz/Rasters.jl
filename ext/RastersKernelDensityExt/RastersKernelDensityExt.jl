module RastersKernelDensityExt

using Extents,
      GeoInterface,
      KernelDensity,
      Rasters

using Rasters: Tables, maybeshiftlocus, Center, lookup, _extent2dims

const GI = GeoInterface

const VectorTuple = Tuple{<:AbstractVector{<:Number},<:AbstractVector{<:Number}}

function Rasters.kerneldensity(geom; geomcolumn=nothing, kw...)
    if Tables.istable(typeof(geom))
        if geomcolumn isa Tuple
            # Use the columns directly and skip _to_vectors
            data = (Tables.getcolumn(geom, geomcolumn[1]), Tables.getcolumn(geom, geomcolumn[2]))
            return _kerneldensity(data; kw...)
        elseif isnothing(geomcolumn)
            geomcolumn = GI.geometrycolumns(geom)[1]
        end
        geom = GI.getcolumn(geom, geomcolumn)
    end
    return _kerneldensity(_to_vectors(geom); kw...)
end
# Also accept the regular KernelDensity.jl data formats
Rasters.kerneldensity(data::VectorTuple; kw...) = _kerneldensity(data; kw...)
Rasters.kerneldensity(data::AbstractMatrix{<:Number}; kw...) =
    _kerneldensity((data[:, 1], data[:, 2]); kw...)

function _kerneldensity(data::VectorTuple;
    to=nothing, size=nothing, res=nothing, crs=nothing, kw...
)
    # Get the extent from `data` if there is no `to` keyword
    if isnothing(to)
        x, y = map(extrema, data)
        to = Extents.Extent(; X=x, Y=y)
    end
    # Convert input keywords to Dimensions
    ds = _extent2dims(to; size, res, crs)

    # KernelDensity needs midpoints
    midpoint_ranges = map(parent ∘ lookup, maybeshiftlocus(Center(), ds))

    # Create a bivariate kernel density estimator
    bkde = KernelDensity.kde(data, midpoint_ranges; kw...)

    # Return a raster
    return Raster(bkde.density, ds)
end

_to_vectors(geoms) = _to_vectors(GI.trait(geoms), geoms)
function _to_vectors(::Nothing, geoms::AbstractArray)
    if all(GI.geomtrait(g) isa PointTrait for g in geoms)
        return _fill_vectors(geoms, length(geoms))
    elseif all(GI.geomtrait(g) isa MultiPointTrait for g in geoms)
        n = sum(GI.npoint(g) for g in geoms)
        points = Iterators.flatten(GI.getpoint(g) for g in geoms)
        return _fill_vectors(points, n)
    end
end
function _to_vectors(::GI.MultiPointTrait, geom)
    n = GI.npoint(geom)
    points = Iterators.flatten(GI.getpoint(g) for g in geoms)
    return _fill_vectors(points, n)
end
_to_vectors(::PointTrait, geom) = [GI.x(geom), GI.y(geom)]
function _to_vectors(trait, geom)
    throw(ArgumentError("$trait not supported: only `MultiPointTrait` `PointTrait` and `MultiPointTrait` are supported"))
end

function _fill_vectors(points, n)
    vec1 = Vector{Float64}(undef, n)
    vec2 = Vector{Float64}(undef, n)
    for (i, p) in enumerate(points)
        vec1[i], vec2[i] = GI.x(p), GI.y(p)
    end
    return (vec1, vec2)
end

end
