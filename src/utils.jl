filter_ext(path, ext::AbstractString) = filter(fn -> splitext(fn)[2] == ext, readdir(path))
filter_ext(path, exts::Union{Tuple,AbstractArray}) =
    filter(fn -> splitext(fn)[2] in exts, readdir(path))
filter_ext(path, ext::Nothing) = readdir(path)

cleankeys(name) = (_cleankey(name),)
function cleankeys(keys::Union{NamedTuple,Tuple,AbstractArray})
    Tuple(map(_cleankey, keys, ntuple(i -> i, length(keys))))
end

function _cleankey(name::Union{Symbol,AbstractString,Name,NoName}, i=1)
    if name in (NoName(), Symbol(""), Name(Symbol("")))
        Symbol("layer$i")
    else
        Symbol(name)
    end
end

noindex_to_sampled(A) = rebuild(A; dims=noindex_to_sampled(dims(A)))
function noindex_to_sampled(dims::DimTuple)
    map(dims) do d
        lookup(d) isa NoLookup ? set(d, Sampled) : d
    end
end

function maybe_typemin_as_missingval(filename::String, A::AbstractRaster{T}) where T
    if ismissing(missingval(A))
        newmissingval = typemin(Missings.nonmissingtype(T))
        base, ext = splitext(filename)
        A1 = replace_missing(A, newmissingval)
        if missing isa eltype(A1)
            A1 = replace_missing(A, missing)
        end
        @warn "`missing` cant be written to $ext, typemin for `$(eltype(A1))` of `$newmissingval` used instead"
        return A1
    elseif missing isa eltype(A)
        A1 = replace_missing(A, missingval)
    else
        return A
    end
end

# Create a standardisted Metadata object of source T, containing a `Dict{String,Any}`
_metadatadict(T::Type, p1::Pair, pairs::Pair...) = _metadatadict(T, (p1, pairs...))
_metadatadict(::Type{T}) where T = Metadata{T}(Dict{String,Any}())
function _metadatadict(::Type{T}, pairs) where T
    dict = Dict{String,Any}()
    for (k, v) in pairs
        dict[String(k)] = v
    end
    return Metadata{T}(dict)
end

# We often need to convert the locus and the lookup in the same step,
# as doing it in the wrong order can give errors.
# function convert_locus_lookup(M1::Type{<:LookupArray}, L1::Type{<:Locus}, dim::Dimension)
#     _convert(S1, L1, sampling(dim), locus(dim), span(dim), dim)
# end

# _convert(::Type{M1}, ::Type{L1}, lookup::M2, l2::L2, span, dim) where {M1,M2<:M1,L1,L2<:L1} = dim
# _convert(::Type{M1}, ::Type{L1}, lookup::M2, l2::L2, span, dim) where {M1,M2<:M1,L1,L2} = shiftlocus(L1(), dim)
# _convert(::Type{M1}, ::Type{L1}, lookup::M2, l2::L2, span, dim) where {M1,M2<:M1,L1,L2} =
#     _convert_by_locus(M1, L1, lookup, l2, span, dim)

# _convert_by_locus(M1, ::Type{Center}, lookup, l2::Union{Start,End}, span, dim) =
#     _convert_by_lookup(M1, dim)
# _convert_by_locus(M1, L1::Type{Union{Start,End}}, lookup, l2::Center, span, dim) =
#     _convert_by_lookup(M1, L1, dim)
# _convert_by_locus(M1, L1::Type{Start}, lookup, l2::End, span, dim) =
#     convertlookup(M1, shiftlocus(L1, dim))
# _convert_by_locus(M1, L1::Type{End}, lookup, l2::Start, span, dim) =
#     convertlookup(M1, shiftlocus(L1, dim))
# _convert_by_locus(M1, ::Type{L1}, lookup, l2::L2, span, dim) where {L1,L2<:L1} =
#     convertlookup(M1, dim)

# # Projected will have an accurate center point on an equal scale, but Mapped may not.
# # So we always shift the locus while in Projected lookup to avoid errors.
# _convert_by_lookup(::Type{Mapped}, dim) = convertlookup(Mapped, shiftlocus(Center(), dim))
# _convert_by_lookup(::Type{Projected}, dim) = shiftlocus(Center(), convertlookup(Projected, dim))


_unwrap(::Val{X}) where X = X
_unwrap(x) = x

_missingval_or_missing(x) = missingval(x) isa Nothing ? missing : missingval(x)

maybe_eps(dims::DimTuple) = map(maybe_eps, dims)
maybe_eps(dim::Dimension) = maybe_eps(eltype(dim))
maybe_eps(::Type) = nothing
maybe_eps(T::Type{<:AbstractFloat}) = _default_atol(T)

_writeable_missing(filename::Nothing, T) = missing
_writeable_missing(filename::AbstractString, T) = _writeable_missing(T)
function _writeable_missing(T)
    missingval = typemin(Missings.nonmissingtype(T))
    @info "`missingval` set to typemin of $missingval"
    return missingval
end

# Map filename suffix over a stack
function mapargs(f, st::AbstractRasterStack, args...)
    layers = map(values(st), args...) do A, mappedargs...
        f(A, mappedargs...)
    end
    return DD.rebuild_from_arrays(st, Tuple(layers))
end

_without_mapped_crs(f, x) = _without_mapped_crs(f, x, mappedcrs(x))
_without_mapped_crs(f, x, ::Nothing) = f(x)
function _without_mapped_crs(f, dims::DimTuple, mappedcrs::GeoFormat)
    dims1 = setmappedcrs(dims, nothing)
    x = f(dims1)
    if x isa DimTuple
        x = setmappedcrs(x, mappedcrs)
    end
    return x
end
function _without_mapped_crs(f, A::AbstractRaster, mappedcrs::GeoFormat)
    A = setmappedcrs(A, nothing)
    x = f(A)
    if x isa AbstractRaster
        x = setmappedcrs(x, mappedcrs)
    end
    return x
end
function _without_mapped_crs(f, st::AbstractRasterStack, mappedcrs::GeoFormat)
    st1 = map(A -> setmappedcrs(A, nothing), st)
    x = f(st1)
    if x isa AbstractRasterStack
        x = map(A -> setmappedcrs(A, mappedcrs(st)), x)
    end
    return x
end

function _extent2dims(to::Extents.Extent{K};
    size=nothing, res=nothing, crs=nothing, kw...
) where K
    emptydims = map(key2dim, K)
    if isnothing(size)
        isnothing(res) && throw(ArgumentError("Pass either `size` or `res` keywords or a `Tuple` of `Dimension`s for `to`."))
        if res isa Real
            res = ntuple(_ -> res, length(K))
        end
        ranges = map(values(to), res) do bounds, r
            start, outer = bounds
            length = ceil(Int, (outer - start) / r)
            step = (outer - start) / length
            range(; start, step, length)
        end
    else
        isnothing(res) || throw(ArgumentError("Both `size` and `res` keywords are passed, but only one can be used"))
        if size isa Int
            size = ntuple(_ -> size, length(K))
        end
        ranges = map(values(to), size) do bounds, length
            start, outer = bounds
            step = (outer - start) / length
            range(; start, step, length)
        end
    end
    lookups = map(ranges) do range
        Projected(range;
            order=ForwardOrdered(),
            sampling=Intervals(Start()),
            span=Regular(step(range)),
            crs,
        )
    end
    d = map(rebuild, emptydims, lookups)
    return d
end


# Like `create` but without disk writes, mostly for Bool/Union{Missing,Boo},
# and uses `similar` where possible
# TODO merge this with `create` somehow
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
    _init_bools(to, dims, T, data; kw...)
end
function _init_bools(to, dims::DimTuple, T::Type, data; combine=true, kw...)
    if combine
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

# When `to` is a Raster we can try to use the same parent array type
function _alloc_bools(to::AbstractRaster, dims::DimTuple, ::Type{T}; missingval=nothing, kw...) where T
    dims = commondims(dims, DEFAULT_POINT_ORDER)
    # TODO: improve this so that only e.g. CuArray uses `similar`
    # This is a little annoying to lock down for all wrapper types,
    # maybe ArrayInterface has tools for this.
    data = if T === Bool && parent(to) isa Union{Array,DA.AbstractDiskArray} 
        # falses(dims) # Use a BitArray
        fill!(similar(to, T, dims), missingval) # Fill some other array type
    else
        fill!(similar(to, T, dims), missingval) # Fill some other array type
    end
    return Raster(data, dims; missingval)
end
# Otherwise just use an Array or BitArray
function _alloc_bools(to, dims::DimTuple, ::Type{T}; missingval, kw...) where T
    data = if T === Bool
        # falses(dims) # Use a BitArray
        fill!(Raster{T}(undef, dims), missingval) # Use an `Array`
    else
        fill!(Raster{T}(undef, dims), missingval) # Use an `Array`
    end
    return Raster(data, dims; missingval)
end

function _as_intervals(ds::Tuple)
    # Rasterization onl makes sense on Intervals
    interval_dims = map(dims(ds, DEFAULT_POINT_ORDER)) do d
        set(d, Sampled(; sampling=Intervals()))
    end
    return setdims(ds, interval_dims)
end


_warn_disk() = @warn "Disk-based objects may be very slow here. User `read` first."

_filenotfound_error(filename) = throw(ArgumentError("file \"$filename\" not found"))

