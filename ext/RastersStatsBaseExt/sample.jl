Rasters.sample(x::RA.RasterStackOrArray, n::Integer; kw...) = Rasters.sample(Random.GLOBAL_RNG, x, n; kw...)
@inline function Rasters.sample(
    rng::Random.AbstractRNG, x::RA.RasterStackOrArray, n::Integer; 
    geometry=true,
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
 #   x = x isa RA.AbstractRasterStack ? x[na] : x
    _sample(rng, x, n;
        dims=DD.dims(x, RA.DEFAULT_POINT_ORDER),
        names=NamedTuple{na}(na),
        # These keywords are converted to _True/_False for type stability later on
        # The @inline above helps constant propagation of the Bools
        geometry=_booltype(geometry), 
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
    dims, names::NamedTuple{K}, geometry, index, skipmissing, weights, replace, ordered, weightstype
) where K
    indices = sample_indices(rng, x, n, skipmissing, weights, replace, ordered, weightstype)
    tuplepoint = map(first, dims)
    T = RA._rowtype(x, tuplepoint; geometry, index, skipmissing, skipinvalid = _True(), names)
    T2 = RA._rowtype(x, tuplepoint; geometry, index, skipmissing = _False(), skipinvalid = _True(), names)
    points = DimPoints(dims)
    rows = Vector{T}(undef, n)
    for (i, I) in enumerate(indices)
        rows[i] = _getindex(T, x, points, I, names)
    end
    rows = T === T2 ? rows : T.(rows)
    return rows
end

_getindex(::Type{T}, x::AbstractRaster{U}, points, idx, ::NamedTuple{K}) where {T, U, K} = 
    RA._maybe_add_fields(T, NamedTuple{K, Tuple{U}}((x[idx],)), points[idx], idx)
_getindex(::Type{T}, x::AbstractRasterStack{<:Any, K}, points, idx, ::NamedTuple{L}) where {T, K, L} = 
    RA._maybe_add_fields(T, K(x[idx])[L], points[idx], idx)

function sample_indices(rng, x, n, skipmissing::_False, weights::Nothing, replace, ordered, weightstype)
    StatsBase.sample(rng, CartesianIndices(x), n; replace, ordered)
end
function sample_indices(rng, x, n, skipmissing::_True, weights::Nothing, replace, ordered, weightstype)
    wts = weightstype(vec(boolmask(x)))
    StatsBase.sample(rng, CartesianIndices(x), wts, n; replace, ordered)
end
function sample_indices(rng, x, n, skipmissing::_False, weights::AbstractDimArray, replace, ordered, weightstype)
    wts = if dims(weights) == dims(x)
        weights
    else
        @d ones(eltype(weights), dims(x)) .* weights
    end |> vec |> weightstype
    StatsBase.sample(rng, CartesianIndices(x), wts, n; replace, ordered)
end
function sample_indices(rng, x, n, skipmissing::_True, weights::AbstractDimArray, replace, ordered, weightstype)
    wts = weightstype(vec(@d boolmask(x) .* weights))
    StatsBase.sample(rng, CartesianIndices(x), wts, n; replace, ordered)
end
