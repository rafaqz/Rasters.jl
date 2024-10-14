Rasters.sample(x::RA.RasterStackOrArray, n::Integer; kw...) = Rasters.sample(Random.GLOBAL_RNG, x, n; kw...)
@inline function Rasters.sample(
    rng::Random.AbstractRNG, x::RA.RasterStackOrArray, n::Integer; 
    geometry=(X,Y),
    index = false, 
    names=RA._names(x), 
    name=names, 
    skipmissing=false,
    replace=true, 
    ordered=false, 
    weights=nothing, 
    weightstype::Type{<:StatsBase.AbstractWeights}=StatsBase.Weights
)
    na = DD._astuple(name)
    geometry, geometrytype, dims = _geometrytype(x, geometry)
 #   x = x isa RA.AbstractRasterStack ? x[na] : x
    _sample(rng, x, n;
        dims=DD.dims(x, RA.DEFAULT_POINT_ORDER),
        names=NamedTuple{na}(na),
        # These keywords are converted to _True/_False for type stability later on
        # The @inline above helps constant propagation of the Bools
        geometry=geometry,
        geometrytype=geometrytype,  
        index=_booltype(index), 
        skipmissing=_booltype(skipmissing), 
        weights,
        replace, 
        ordered,
        weightstype
    )

end
function _sample(
    rng, x, n; 
    dims, names::NamedTuple{K}, geometry, geometrytype, index, skipmissing, weights, replace, ordered, weightstype
) where K
    indices = sample_indices(rng, x, n, skipmissing, weights, replace, ordered, weightstype)
    points = DimPoints(dims)
    T = RA._rowtype(x, geometrytype; geometry, index, skipmissing, skipinvalid = _True(), names)
    x2 = x isa AbstractRasterStack ? x[K] : RasterStack(NamedTuple{K}((x,)))
    _getindices(T, x2, points, indices)
end
function _getindices(::Type{T}, x, points, indices) where T
    rows = Vector{T}(undef, size(indices))
    for (i, I) in enumerate(indices)
        rows[i] = _getindex(T, x, points, I)
    end
    return rows
end

_getindex(::Type{T}, x::AbstractRasterStack{<:Any, NT}, points, idx) where {T, NT} = 
    RA._maybe_add_fields(T, NT(x[RA.commondims(idx, x)]), points[RA.commondims(idx, points)], val(idx))

function sample_indices(rng, x, n, skipmissing::_False, weights::Nothing, replace, ordered, weightstype)
    StatsBase.sample(rng, RA.DimIndices(x), n; replace, ordered)
end
function sample_indices(rng, x, n, skipmissing::_True, weights::Nothing, replace, ordered, weightstype)
    wts = weightstype(vec(boolmask(x)))
    StatsBase.sample(rng, RA.DimIndices(x), wts, n; replace, ordered)
end
function sample_indices(rng, x, n, skipmissing::_False, weights::AbstractDimArray, replace, ordered, weightstype)
    wts = if dims(weights) == dims(x)
        weights
    else
        @d ones(eltype(weights), dims(x)) .* weights
    end |> vec |> weightstype
    StatsBase.sample(rng, RA.DimIndices(x), wts, n; replace, ordered)
end
function sample_indices(rng, x, n, skipmissing::_True, weights::AbstractDimArray, replace, ordered, weightstype)
    wts = weightstype(vec(@d boolmask(x) .* weights))
    StatsBase.sample(rng, RA.DimIndices(x), wts, n; replace, ordered)
end

function _geometrytype(x, geometry::Bool)
    if geometry
        error("Specify a geometry type by setting `geometry` to a Tuple or NamedTuple of Dimensions. E.g. `geometry = (X, Y)`")
    else
        return _False(), Nothing, dims(x)
    end
end

function _geometrytype(x, geometry::Tuple)
    dims = DD.commondims(DD.dims(x), geometry)
    return _True(), Tuple{map(eltype, dims)...}, dims
end
function _geometrytype(x, geometry::NamedTuple{K}) where K
    dims = DD.commondims(DD.dims(x), values(geometry))
    return _True(), NamedTuple{K, Tuple{map(eltype, dims)...}}, dims
end