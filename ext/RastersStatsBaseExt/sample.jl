Rasters.sample(x::RA.RasterStackOrArray, n::Integer; kw...) = Rasters.sample(Random.GLOBAL_RNG, x, n; kw...)
@inline function Rasters.sample(
    rng::Random.AbstractRNG, x::RA.RasterStackOrArray, n::Integer; 
    geometry = true, index = false, names=RA._names(x), name=names, skipmissing = false,
    replace = true, ordered = false, weights = nothing, weightstype::Type{<:StatsBase.AbstractWeights} = StatsBase.Weights
)
    na = DD._astuple(name)
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
    T = _srowtype(x, tuplepoint; geometry, index, skipmissing, names)
    rows = Vector{T}(undef, n)
    points = DimPoints(dims)
    
    for i in 1:n
        idx = indices[i]
        props = if x isa Raster 
            NamedTuple{K,Tuple{eltype(x)}}((x[idx],))
        else
            NamedTuple(x[idx])[K]
        end
        point = geometry isa _True ? points[idx] : nothing
        rows[i] = RA._maybe_add_fields(T, props, point, idx)
    end

    return rows
end

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

# Determine the row type, making use of some of extract machinery
_srowtype(x, g; kw...) = _srowtype(x, typeof(g); kw...)
function _srowtype(x, g::Type; geometry, index, skipmissing, names, kw...)
    keys = RA._rowkeys(geometry, index, names)
    types = _srowtypes(x, g, geometry, index, skipmissing, names)
    NamedTuple{keys,types}
end
function _srowtypes(
    x, ::Type{G}, geometry::_True, index::_True, skipmissing::_False, names::NamedTuple{Names}
) where {G,Names}
    Tuple{G,Tuple{Int,Int},_nametypes(x, names)...}
end
function _srowtypes(
    x, ::Type{G}, geometry::_False, index::_True, skipmissing::_False, names::NamedTuple{Names}
) where {G,Names}
    Tuple{Tuple{Int,Int},_nametypes(x, names)...}
end
function _srowtypes(
    x, ::Type{G}, geometry::_True, index::_False, skipmissing::_False, names::NamedTuple{Names}
) where {G,Names}
    Tuple{G,_nametypes(x, names)...}
end
function _srowtypes(
    x, ::Type{G}, geometry::_False, index::_False, skipmissing::_False, names::NamedTuple{Names}
) where {G,Names}
    Tuple{_nametypes(x, names)...}
end
# fallback
_srowtypes(x, T, geometry, index, skipmissing::_True, names) = 
    RA._rowtypes(x, T, geometry, index, skipmissing, names)
# adapted from extract code
@inline _nametypes(::Raster{T}, ::NamedTuple{Names}) where {T,Names} = (T,)
function _nametypes(::RasterStack{<:Any,T}, ::NamedTuple{PropNames}) where 
        {T<:NamedTuple{StackNames,Types},PropNames} where {StackNames,Types}
    nt = NamedTuple{StackNames}(Types.parameters)
    return values(nt[PropNames])
end

