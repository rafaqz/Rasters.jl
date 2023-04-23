filter_ext(path, ext::AbstractString) =
    filter(fn -> splitext(fn)[2] == ext, readdir(path; join=true))
filter_ext(path, exts::Union{Tuple,AbstractArray}) =
    filter(fn -> splitext(fn)[2] in exts, readdir(path; join=true))
filter_ext(path, ext::Nothing) = readdir(path; join=true)

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

function _maybe_use_type_missingval(filename::String, A::AbstractRaster{T}) where T
    if ismissing(missingval(A))
        newmissingval = _type_missingval(Missings.nonmissingtype(T))
        base, ext = splitext(filename)
        A1 = replace_missing(A, newmissingval)
        @warn "`missing` cant be written to $ext, missinval for `$(eltype(A1))` of `$newmissingval` used instead"
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
    missingval = _type_missingval(Missings.nonmissingtype(T))
    @info "`missingval` set to $missingval"
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

function _extent2dims(to; size=nothing, res=nothing, crs=nothing, kw...) 
    _extent2dims(to, size, res, crs)
end
function _extent2dims(to::Extents.Extent, size::Nothing, res::Nothing, crs)
    isnothing(res) && throw(ArgumentError("Pass either `size` or `res` keywords or a `Tuple` of `Dimension`s for `to`."))
end
function _extent2dims(to::Extents.Extent, size, res, crs)
    isnothing(res) || _size_and_res_error()
end
function _extent2dims(to::Extents.Extent{K}, size::Nothing, res::Real, crs) where K
    tuple_res = ntuple(_ -> res, length(K))
    _extent2dims(to, size, tuple_res, crs)
end
function _extent2dims(to::Extents.Extent, size::Nothing, res, crs)
    ranges = map(values(to), res) do bounds, r
        start, outer = bounds
        length = ceil(Int, (outer - start) / r)
        step = (outer - start) / length
        range(; start, step, length)
    end
    return _extent2dims(to, ranges, crs)
end
function _extent2dims(to::Extents.Extent{K}, size, res::Nothing, crs) where K
    if size isa Int
        size = ntuple(_ -> size, length(K))
    end
    ranges = map(values(to), size) do bounds, length
        start, outer = bounds
        step = (outer - start) / length
        range(; start, step, length)
    end
    return _extent2dims(to, ranges, crs)
end
function _extent2dims(to::Extents.Extent{K}, ranges, crs) where K
    emptydims = map(key2dim, K)
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

function _as_intervals(ds::Tuple)
    # Rasterization only makes sense on Sampled Intervals
    interval_dims = map(dims(ds, DEFAULT_POINT_ORDER)) do d
        l = parent(d)
        rebuild(d, rebuild(l; sampling=Intervals(locus(l))))
    end
    return setdims(ds, interval_dims)
end

_geomindices(geoms) = GI.isfeaturecollection(geoms) ? (1:GI.nfeature(geoms)) : eachindex(geoms)
_getgeom(geoms, i::Integer) = GI.isfeaturecollection(geoms) ? GI.getfeature(geoms, i) : geoms[i]


_warn_disk() = @warn "Disk-based objects may be very slow here. User `read` first."

_filenotfound_error(filename) = throw(ArgumentError("file \"$filename\" not found"))

_progress(args...; kw...) = ProgressMeter.Progress(args...; color=:blue, barlen=50, kw...)

# Function barrier for splatted vector broadcast
@noinline _do_broadcast!(f, x, args...) = broadcast!(f, x, args...)

_size_and_res_error() = throw(ArgumentError("Both `size` and `res` keywords are passed, but only one can be used"))

_type_missingval(::Type{T}) where T = typemin(T)
_type_missingval(::Type{T}) where T<:Unsigned = typemax(T) 

# Modified from IsURL.jl, many thanks to @zlatanvasovic
const WINDOWSREGEX = r"^[a-zA-Z]:[\\]"
const URLREGEX = r"^[a-zA-Z][a-zA-Z\d+\-.]*:"

_isurl(str::AbstractString) = !occursin(WINDOWSREGEX, str) && occursin(URLREGEX, str)
