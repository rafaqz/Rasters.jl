
# Like `create` but without disk writes, mostly for Bool or Union{Missing,Bool},
# and uses `similar` where possible
# TODO merge this with `create` somehow
_init_bools(to; kw...) = _init_bools(to, BitArray; kw...)
_init_bools(to, T::Type; kw...) = _init_bools(to, T, nothing; kw...)
_init_bools(to::AbstractRasterSeries, T::Type, data; kw...) = _init_bools(first(to), T, data; kw...)
_init_bools(to::AbstractRasterStack, T::Type, data; kw...) = _init_bools(first(to), T, data; kw...)
_init_bools(to::AbstractRaster, T::Type, data; kw...) = _init_bools(to, dims(to), T, data; kw...)
_init_bools(to::DimTuple, T::Type, data; kw...) = _init_bools(to, to, T, data; kw...)
function _init_bools(to::Nothing, T::Type, data; 
    geometrycolumn=nothing, 
    collapse=nokw, 
    res=nokw,
    size=nokw,
    kw...
)
    # Get the extent of the geometries
    ext = _extent(data; geometrycolumn)
    isnothing(ext) && throw(ArgumentError("no recognised dimensions, extent or geometry"))
    return _init_bools(ext, T, data; collapse, res, size)
end
function _init_bools(to::Extents.Extent, T::Type, data;
    collapse=nokw, size=nokw, res=nokw, sampling=nokw, kw...
)
    # Convert the extent to dims (there must be `res` or `size` in `kw`)
    ext = _extent2dims(to; size, res, sampling, kw...)
    _init_bools(to, ext, T, data; collapse, kw...)
end
function _init_bools(to, dims::DimTuple, T::Type, data; 
    collapse::Union{Bool,Nothing,NoKW}=nokw, kw...
)
    if isnothing(data) || isnokwornothing(collapse) || collapse
        _alloc_bools(to, dims, T; kw...)
    else
        n = if Base.IteratorSize(data) isa Base.HasShape
            length(data)
        else
            count(_ -> true, data)
        end
        geomdim = Dim{:geometry}(1:n)
        _alloc_bools(to, (dims..., geomdim), T; kw...)
    end
end

function _alloc_bools(to, dims::DimTuple, ::Type{BitArray}; missingval::Bool=false, metadata=NoMetadata(), kw...)
    # Use a BitArray
    vals = missingval == false ? falses(size(dims)) : trues(size(dims))
    return Raster(vals, dims; missingval, metadata)
end
function _alloc_bools(to, dims::DimTuple, ::Type{<:Array{T}}; missingval=false, metadata=NoMetadata(), kw...) where T
    # Use an Array
    data = fill!(Raster{T}(undef, dims), missingval) 
    return rebuild(data; missingval, metadata)
end

function _prepare_for_burning(B, locus=Center())
    B1 = _forward_ordered(B)
    start_dims = map(dims(B1, DEFAULT_POINT_ORDER)) do d
        # Shift lookup values to center of pixels
        d = DD.maybeshiftlocus(locus, d)
        _lookup_as_array(d)
    end
    return setdims(B1, start_dims)
end

# Convert to Array if its not one already
_lookup_as_array(d::Dimension) = parent(lookup(d)) isa Array ? d : modify(Array, d) 

function _forward_ordered(B)
    reduce(dims(B); init=B) do A, d
        if DD.order(d) isa ReverseOrdered
            A = view(A, rebuild(d, lastindex(d):-1:firstindex(d)))
            set(A, d => reverse(d))
        else
            A
        end
    end
end

_nthreads() = Threads.nthreads()
