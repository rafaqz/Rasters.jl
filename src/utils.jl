
# File paths, urls and strings

filter_ext(path, ext::AbstractString) =
    filter(fn -> splitext(fn)[2] == ext, readdir(path; join=true))
filter_ext(path, exts::Union{Tuple,AbstractArray}) =
    filter(fn -> splitext(fn)[2] in exts, readdir(path; join=true))
filter_ext(path, ext::Nothing) = readdir(path; join=true)

_maybe_add_suffix(filename::Nothing, suffix) = nothing
_maybe_add_suffix(filename::Nothing, suffix::Union{Nothing,NoKW}) = nothing
_maybe_add_suffix(filename, suffix::Union{Nothing,NoKW}) = filename
function _maybe_add_suffix(filename, suffix)
    base, ext = splitext(filename)
    if string(suffix) == ""
        filename
    else
        return string(base, "_", suffix, ext)
    end
end

# Modified from IsURL.jl, many thanks to @zlatanvasovic
const WINDOWSREGEX = r"^[a-zA-Z]:[\\]"
const URLREGEX = r"^[a-zA-Z][a-zA-Z\d+\-.]*:"

_isurl(str::AbstractString) = !occursin(WINDOWSREGEX, str) && occursin(URLREGEX, str)

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

# We often need to convert the locus and the lookup in the same step,
# as doing it in the wrong order can give errors.
# function convert_locus_lookup(M1::Type{<:Lookup}, L1::Type{<:Locus}, dim::Dimension)
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


# Missing values

_missingval_or_missing(x) = _maybe_to_missing(missingval(x))

_maybe_to_missing(::Union{Nothing,NoKW}) = missing
_maybe_to_missing(missingval) = missingval

_writeable_missing(filename::Nothing, T; kw...) = missing
_writeable_missing(filename::AbstractString, T; kw...) = _writeable_missing(T; kw...)
function _writeable_missing(::Type{Missing}; verbose=true)
    missingval = _type_missingval(UInt8)
    verbose && @info "`missingval` set to $missingval on disk"
    return missingval
end
function _writeable_missing(T; verbose=true)
    missingval = _type_missingval(Missings.nonmissingtype(T))
    verbose && @info "`missingval` set to $missingval on disk"
    return missingval
end

_type_missingval(::Type{T}) where T = _type_missingval1(Missings.nonmissingtype(T))
_type_missingval1(::Type{T}) where T = typemin(T)
_type_missingval1(::Type{T}) where T<:Unsigned = typemax(T)

_fix_missingval(::Type, ::Union{NoKW,Nothing}) = nothing
_fix_missingval(::AbstractArray, ::Nothing) = nothing
_fix_missingval(A::AbstractArray, ::NoKW) = _fix_missingval(A, Rasters.missingval(A))
_fix_missingval(::AbstractArray{T}, missingval) where T = _fix_missingval(T, missingval)
function _fix_missingval(::Type{T}, missingval::M) where {T,M}
    T1 = nonmissingtype(T)
    if missingval isa T
        missingval
    elseif hasmethod(convert, Tuple{Type{T1},M}) && isreal(missingval) && 
            missingval <= typemax(T1) && missingval >= typemin(T1)
        if T1 <: Integer && !isinteger(missingval) 
            nothing
        else
            convert(T, missingval)
        end
    else
        nothing
    end
end


# Extents

function _extent2dims(to::Extents.Extent; 
    size=nokw, res=nokw, crs=nokw, mappedcrs=nokw, sampling=nokw,
)
    sampling = _match_to_extent(to, isnokw(sampling) ? Intervals(Start()) : sampling)
    _extent2dims(to, size, res; crs, mappedcrs, sampling)
end
function _extent2dims(to::Extents.Extent, size::Union{Nothing,NoKW}, res::Union{Nothing,NoKW}; kw...)
    throw(ArgumentError("Pass either `size` or `res` keywords or a `Tuple` of `Dimension`s for `to`."))
end
_extent2dims(to::Extents.Extent, size, res; kw...) = _size_and_res_error()
function _extent2dims(to::Extents.Extent, size::Union{Nothing,NoKW}, res; kw...)
    _extent2dims(to, size, _match_to_extent(to, res); kw...)
end
function _extent2dims(to::Extents.Extent, size::Union{Nothing,NoKW}, res::Tuple; 
    sampling::Tuple, kw...
)
    ranges = map(values(to), res, sampling) do (start, stop), step, s
        r = if step >= zero(step)
            range(; start, step, stop)
        else
            range(; start=stop, step, stop=start)
        end
        r = if s isa Intervals
            if locus(s) isa Start
                r[1:end-1]
            elseif locus(s) isa End
                r[2:end]
            else # Center
                r .+ abs(step) / 2
            end
        end
    end
    return _extent2dims(to, ranges; sampling, kw...)
end
_extent2dims(to::Extents.Extent, size, res::Union{Nothing,NoKW}; kw...) =
    _extent2dims(to, _match_to_extent(to, size), res; kw...)
function _extent2dims(to::Extents.Extent, size::Tuple, res::Union{Nothing,NoKW};
    sampling::Tuple, crs, mappedcrs
)
    ranges = map(values(to), size, sampling) do (start, stop), length, s
        if s isa Points
            range(; start, stop, length)
        else
            range(; start, stop, length=length+1)[1:end-1]
        end
    end
    return _extent2dims(to, ranges; sampling, crs, mappedcrs)
end
function _extent2dims(::Extents.Extent{K}, ranges; crs, mappedcrs, sampling::Tuple) where K
    crs = isnokw(crs) ? nothing : crs 
    mappedcrs = isnokw(mappedcrs) ? nothing : mappedcrs 
    emptydims = map(name2dim, K)
    lookups = map(emptydims, ranges, sampling) do d, range, s
        order = Lookups.orderof(range)
        span = Regular(step(range))
        if d isa SpatialDim && !isnothing(crs)
            Projected(range; sampling=s, order, span, crs, mappedcrs)
        else
            Sampled(range; sampling=s, order, span)
        end
    end
    d = map(rebuild, emptydims, lookups)
    return d
end

function _match_to_extent(::Extents.Extent{K}, x) where K
    if x isa DimTuple 
        map(val, dims(x, map(name2dim, K)))
    elseif x isa NamedTuple 
        values(x[K])
    elseif x isa Tuple 
        x
    else
        map(_ -> x, K)
    end
end

# Geometries

# get geometries from what may be a table with a geometrycolumn or an interable of geometries
# if it has no geometry column and does not iterate valid geometries, error informatively
function _get_geometries(data, ::Nothing)
    # if it's a table, get the geometry column
    geoms = if !(data isa AbstractVector{<:GeoInterface.NamedTuplePoint}) && Tables.istable(data)
        geomcol = first(GI.geometrycolumns(data))
        !in(geomcol, Tables.columnnames(Tables.columns(data))) &&
            throw(ArgumentError("Expected geometries in the column `$geomcol`, but no such column found."))
        isnothing(geomcol) && throw(ArgumentError("No default `geometrycolumn` for this type, please specify it manually."))
        Tables.getcolumn(Tables.columns(data), geomcol)
    elseif data isa AbstractVector
        data
    else
        trait = GI.trait(data)
        if GI.trait(data) isa GI.FeatureCollectionTrait
            [GI.geometry(f) for f in GI.getfeature(data)]
        else
            collect(data)
        end
    end
    # check if data iterates valid geometries before returning
    _check_geometries(geoms)
    return geoms
end
function _get_geometries(data, geometrycolumn::Symbol)
    Tables.istable(data) || throw(ArgumentError("`geometrycolumn` was specified, but `data` is not a table."))
    geoms = Tables.getcolumn(Tables.columns(data), geometrycolumn)
    _check_geometries(geoms)
    return geoms
end
function _get_geometries(data, geometrycolumn::NTuple{<:Any, <:Symbol})
    Tables.istable(data) || throw(ArgumentError("`geometrycolumn` was specified, but `data` is not a table."))
    cols = Tables.columns(data)
    geomcols = (Tables.getcolumn(cols, col) for col in geometrycolumn)
    points = map(geomcols...) do (row...)
        for r in row
            ismissing(r) && return missing
        end
        return row
    end
    return points
end
function _check_geometries(geoms)
    for g in geoms
        ismissing(g) || GI.geomtrait(g) !== nothing ||
        throw(ArgumentError("$g is not a valid GeoInterface.jl geometry"))
    end
    return
end
# to distinguish between objects returned by _get_geometries and other objects
struct IterableOfGeometries end


# Chunking

# NoKW means true
@inline function _chunks_to_tuple(template, dimorder, chunks::Bool)
    if chunks == true
        if template isa AbstractArray && DA.haschunks(template) == DA.Chunked()
            # Get chunks from the template
            DA.max_chunksize(DA.eachchunk(template))
        else
            # Use defaults
            _chunks_to_tuple(template, dimorder, (X(512), Y(512)))
        end
    else
        nothing
    end
end
@inline function _chunks_to_tuple(template, dimorder, chunks::NTuple{N,Integer}) where N
    n = length(dimorder)
    if n < N
        throw(ArgumentError("Length $n tuple needed for `chunks`, got $N"))
    elseif n > N
        (chunks..., ntuple(_ -> 1, Val{n-N}())...)
    else
        chunks
    end
end
@inline function _chunks_to_tuple(template, dimorder, chunks::DimTuple)
    size_one_chunk_axes = map(d -> rebuild(d, 1), otherdims(dimorder, chunks))
    alldims = (chunks..., size_one_chunk_axes...)
    int_chunks = map(val, dims(alldims, dimorder))
    if !isnothing(template)
        if !all(map(>=, size(template), int_chunks))
            @warn "Chunks $int_chunks larger than array size $(size(template)). Using defaults."
            return nothing
        end
    end
    return int_chunks
end
@inline _chunks_to_tuple(template, dimorder, chunks::NamedTuple) =
    _chunks_to_tuple(template, dimorder, DD.kw2dims(chunks))
@inline _chunks_to_tuple(template, dimorder, chunks::Nothing) = nothing
@inline _chunks_to_tuple(template, dims, chunks::NoKW) = nothing

_checkregular(A::AbstractRange) = true
function _checkregular(A::AbstractArray)
    step = stepof(A)
    for i in eachindex(A)[2:end]
        if !(A[i] - A[i-1] ≈ step)
            return false
        end
    end
    return true
end

# Constructor helpers

function _raw_check(raw, scaled, maskingval)
    if raw
        scaled isa Bool && scaled && @warn "`scaled=true` set to `false` because of `raw=true`"
        isnokwornothing(maskingval) || @warn "`maskingval=$maskingval` set to `nothing` because of `raw=true`"
        return false, nothing
    else
        scaled = isnokw(scaled) ? true : scaled
        return scaled, maskingval 
    end
end

function _maybe_drop_single_band(x, dropband::Bool, lazy::Bool)
    dropband || return x
    if hasdim(x, Band()) && size(x, Band()) < 2
         if lazy
             return view(x, Band(1)) # TODO fix dropdims in DiskArrays
         else
             return dropdims(x; dims=Band())
         end
    else
         return x
    end
end


# Memory

function _checkobjmem(obj)
    f = bytes -> """
        required memory $(bytes) is greater than system memory $(Sys.free_memory()).
        Use `lazy=true` if you are loading dataset, and only call `read` on a subset after `view`.
        """
    _checkobjmem(f, obj)
end
_checkobjmem(f, obj) = _checkmem(f, _sizeof(obj))

_checkmem(f, bytes::Int) = Sys.free_memory() > bytes || _no_memory_error(f, bytes)

_sizeof(A::AbstractArray{T}) where T = sizeof(T) * prod(size(A))
_sizeof(st::AbstractRasterStack) = sum(_sizeof, layers(st))
_sizeof(s::AbstractRasterSeries) =
    length(s) == 0 ? 0 : _sizeof(first(s)) * prod(size(s))

function _no_memory_error(f, bytes)
    msg = f(bytes) * """
    If you beleive this is not correct, pass the keyword `checkmem=false` or set `Rasters.checkmem!(false)`
    and try again. These options may crash your system if the file is actually larger than memory.
    """
    return error(msg)
end


# Lookups

function _as_intervals(ds::Tuple)
    # Rasterization only makes sense on Sampled Intervals
    interval_dims = map(dims(ds, DEFAULT_POINT_ORDER)) do d
        l = parent(d)
        rebuild(d, rebuild(l; sampling=Intervals(locus(l))))
    end
    return setdims(ds, interval_dims)
end

nolookup_to_sampled(A) = rebuild(A; dims=nolookup_to_sampled(dims(A)))
nolookup_to_sampled(dims::DimTuple) = map(nolookup_to_sampled, dims)
nolookup_to_sampled(d::Dimension) =
    lookup(d) isa NoLookup ? set(d, Sampled(; sampling=Points())) : d


# Metadata

# Create a standardised Metadata object of source T, containing a `Dict{String,Any}`
_metadatadict(s::Source, p1::Pair, pairs::Pair...) =
    _metadatadict(s, (p1, pairs...))
_metadatadict(::S) where S<:Source = Metadata{S}(Dict{String,Any}())
function _metadatadict(::S, pairs) where S<:Source
    dict = Dict{String,Any}()
    for (k, v) in pairs
        dict[String(k)] = v
    end
    return Metadata{S}(dict)
end


# Other

_progress(args...; kw...) = ProgressMeter.Progress(args...; color=:blue, barlen=50, kw...)

# Function barrier for splatted vector broadcast
@noinline _do_broadcast!(f, x, args...) = broadcast!(f, x, args...)

# Run `f` threaded or not, w
function _run(f, range::OrdinalRange, threaded::Bool, progress::Bool, desc::String)
    p = progress ? _progress(length(range); desc) : nothing
    if threaded
        Threads.@threads :static for i in range
            f(i)
            isnothing(p) || ProgressMeter.next!(p)
        end
    else
        for i in range
            f(i)
            isnothing(p) || ProgressMeter.next!(p)
        end
    end
end

_unwrap(::Val{X}) where X = X
_unwrap(x) = x

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


# Warnings and erros

_maybe_warn_replace_missing(replace_missing::NoKW) = nothing
function _maybe_warn_replace_missing(replace_missing)
    @warn "`replace_missing` keyword no longer used. Set `maskingval` to nothing for no replacement, to `missing` to mask `missingval` with `missing`, or any other value"
end

@noinline _warn_disk() = @warn "Disk-based objects may be very slow here. User `read` first."

_filenotfound_error(filename) = throw(ArgumentError("file \"$filename\" not found"))

_size_and_res_error() = throw(ArgumentError("Both `size` and `res` keywords are passed, but only one can be used"))

_no_crs_error() = throw(ArgumentError("The provided object does not have a CRS. Use `setcrs` to set one."))
