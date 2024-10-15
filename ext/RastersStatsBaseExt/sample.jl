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

    _sample(rng, x, n,
        dims,
        NamedTuple{na}(na),
        geometry,
        geometrytype,  
        # These keywords are converted to _True/_False for type stability later on
        _booltype(index), 
        _booltype(skipmissing), 
        weights,
        weightstype,
        replace, 
        ordered
    )

end
function _sample(
    rng, x, n, 
    dims, names::NamedTuple{K}, geometry, ::Type{G}, index, skipmissing, weights, weightstype, replace, ordered, 
) where {K, G}
    indices = sample_indices(rng, x, n, skipmissing, weights, replace, ordered, weightstype)
    T = RA._rowtype(x, G; geometry, index, skipmissing, skipinvalid = _True(), names)
    x2 = x isa AbstractRasterStack ? x[K] : RasterStack(NamedTuple{K}((x,)))
    return _getindices(T, x2, dims, indices)
end

_getindices(::Type{T}, x, dims, indices) where {T} = 
    broadcast(I -> _getindex(T, x, dims, I), indices)

function _getindex(::Type{T}, x::AbstractRasterStack{<:Any, NT}, dims, idx) where {T, NT}
    RA._maybe_add_fields(
        T, 
        NT(x[RA.commondims(idx, x)]), 
        DimPoints(dims)[RA.commondims(idx, dims)], 
        val(idx)
    )
end

function sample_indices(rng, x, n, skipmissing, weights::Nothing, replace, ordered, _)
    if istrue(skipmissing)
        wts = StatsBase.Weights(vec(boolmask(x)))
        StatsBase.sample(rng, RA.DimIndices(x), wts, n; replace, ordered)
    else
        StatsBase.sample(rng, RA.DimIndices(x), n; replace, ordered)
    end
end
function sample_indices(rng, x, n, skipmissing, weights::AbstractDimArray, replace, ordered, ::Type{W}) where W
    wts = if istrue(skipmissing) 
        @d boolmask(x) .* weights
    elseif dims(weights) == dims(x)
        weights
    else
        @d ones(eltype(weights), dims(x)) .* weights
    end |> vec |> W
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