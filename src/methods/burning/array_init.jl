
# Like `create` but without disk writes, mostly for Bool/Union{Missing,Boo},
# and uses `similar` where possible
# TODO merge this with `create` somehow
_init_bools(to; kw...) = _init_bools(to, Bool; kw...)
_init_bools(to, T::Type; kw...) = _init_bools(to, T, nothing; kw...)
_init_bools(to::AbstractRasterSeries, T::Type, data; kw...) = _init_bools(first(to), T, data; kw...)
_init_bools(to::AbstractRasterStack, T::Type, data; kw...) = _init_bools(first(to), T, data; kw...)
_init_bools(to::AbstractRaster, T::Type, data; kw...) = _init_bools(to, dims(to), T, data; kw...)
_init_bools(to::Extents.Extent, T::Type, data; kw...) = _init_bools(to, _extent2dims(to; kw...), T, data; kw...)
_init_bools(to::DimTuple, T::Type, data; kw...) = _init_bools(to, to, T, data; kw...)
function _init_bools(to::Nothing, T::Type, data; kw...)
    # Get the extent of the geometries
    ext = _extent(data)
    isnothing(ext) && throw(ArgumentError("no recognised dimensions, extent or geometry"))
    # Convert the extent to dims (there must be `res` or `size` in `kw`)
    dims = _extent2dims(ext; kw...)
    return _init_bools(to, dims, T, data; kw...)
end
function _init_bools(to, dims::DimTuple, T::Type, data; collapse::Union{Bool,Nothing}=nothing, kw...)
    if isnothing(data) || isnothing(collapse) || collapse
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

function _alloc_bools(to, dims::DimTuple, ::Type{Bool}; missingval=false, metadata=NoMetadata(), kw...)
    if length(dims) > 2
        # Use a BitArray
        return Raster(falses(size(dims)), dims; missingval, metadata) # Use a BitArray
    else
        return Raster(zeros(Bool, size(dims)), dims; missingval, metadata) # Use a BitArray
    end
end
function _alloc_bools(to, dims::DimTuple, ::Type{T}; missingval=false, metadata=NoMetadata(), kw...) where T
    # Use an `Array`
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
