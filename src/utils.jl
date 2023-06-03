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
function _extent2dims(to::Extents.Extent{K}, size::Nothing, res, crs) where K
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

_geomindices(geoms) = _geomindices(GI.trait(geoms), geoms) 
_geomindices(::Nothing, geoms) = eachindex(geoms)
_geomindices(::GI.FeatureCollectionTrait, geoms) = 1:GI.nfeature(geoms)
_geomindices(::GI.FeatureTrait, geoms) = _geomindices(GI.geometry(geoms))
_geomindices(::GI.AbstractGeometryTrait, geoms) = 1:GI.ngeom(geoms)

_getgeom(geoms, i::Integer) = _getgeom(GI.trait(geoms), geoms, i)
_getgeom(::GI.FeatureCollectionTrait, geoms, i::Integer) = GI.geometry(GI.getfeature(geoms, i))
_getgeom(::GI.FeatureTrait, geoms, i::Integer) = GI.getgeom(GI.geometry(geoms), i)
_getgeom(::GI.AbstractGeometryTrait, geoms, i::Integer) = GI.getgeom(geom, i)
_getgeom(::GI.PointTrait, geom, i::Integer) = error("PointTrait should not be reached")
_getgeom(::Nothing, geoms, i::Integer) = geoms[i] # Otherwise we can probably just index?

_getgeom(geoms) = _getgeom(GI.trait(geoms), geoms)
_getgeom(::GI.FeatureCollectionTrait, geoms) = (GI.geometry(f) for f in GI.getfeature(geoms))
_getgeom(::GI.FeatureTrait, geoms) = GI.getgeom(GI.geometry(geoms))
_getgeom(::GI.AbstractGeometryTrait, geoms) = GI.getgeom(geoms)
_getgeom(::GI.PointTrait, geom) = error("PointTrait should not be reached")
_getgeom(::Nothing, geoms) = geoms


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

_checkbounds(A::AbstractRasterStack, I...) = checkbounds(Bool, first(A), I...)
_checkbounds(A::AbstractRaster, I...) = checkbounds(Bool, A, I...)

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


const XYExtent = Extents.Extent{(:X,:Y),Tuple{Tuple{Float64,Float64},Tuple{Float64,Float64}}}

# Get the bounds of a geometry
_extent(geom; kw...)::XYExtent = _extent(GI.trait(geom), geom; kw...)
function _extent(::Nothing, data::AbstractVector; kw...)::XYExtent
    g1 = first(data)
    if GI.trait(g1) isa GI.PointTrait 
        xs = extrema(p -> GI.x(p), data)
        ys = extrema(p -> GI.y(p), data)
        return _float64_xy_extent(Extents.Extent(X=xs, Y=ys))
    else
        ext = reduce(data; init=_extent(first(data))) do ext, geom
            Extents.union(ext, _extent(geom))
        end
        return _float64_xy_extent(ext)
    end
end
_extent(::Nothing, data::RasterStackOrArray; kw...)::XYExtent = _float64_xy_extent(Extents.extent(data))
function _extent(::Nothing, data::T; geometrycolumn=nothing)::XYExtent where T
    if Tables.istable(T)
        singlecolumn = isnothing(geometrycolumn) ? first(GI.geometrycolumns(data)) : geometrycolumn
        cols = Tables.columns(data)
        if singlecolumn isa Symbol && singlecolumn in Tables.columnnames(cols)
            # Table of geometries
            geoms = Tables.getcolumn(data, singlecolumn)
            return _extent(nothing, geoms)
        else
            multicolumn = isnothing(geometrycolumn) ? DEFAULT_POINT_ORDER : geometrycolumn 
            # TODO: test this branch
            # Table of points with dimension columns
            bounds = reduce(multicolumn; init=(;)) do acc, key
                if key in Tables.columnnames(cols)
                    merge(acc, (; key=extrema(cols[key])))
                else
                    acc
                end
            end
            return _float64_xy_extent(Extensts.Extent(bounds))
        end
    else
        ext = Extents.extent(data)
        ext isa Extents.Extent || throw(ArgumentError("object returns `nothing` from `Extents.extent`."))
        return _float64_xy_extent(ext)
    end
end
function _extent(::GI.AbstractPointTrait, point; kw...)::XYExtent
    x, y = Float64(GI.x(point)), Float64(GI.y(point))
    Extents.Extent(X=(x, x), Y=(y, y))
end
function _extent(::GI.AbstractGeometryTrait, geom; kw...)::XYExtent
    geomextent = GI.extent(geom; fallback=false)
    if isnothing(geomextent)
        points = GI.getpoint(geom)
        xbounds = extrema(GI.x(p) for p in points)
        ybounds = extrema(GI.y(p) for p in points)
        return _float64_xy_extent(Extents.Extent(X=xbounds, Y=ybounds))
    else
        return _float64_xy_extent(geomextent)
    end
end
_extent(::GI.AbstractFeatureTrait, feature; kw...)::XYExtent = _extent(GI.geometry(feature))
function _extent(::GI.AbstractFeatureCollectionTrait, features; kw...)::XYExtent
    features = GI.getfeature(features)
    init = _float64_xy_extent(_extent(first(features)))
    ext = reduce(features; init) do acc, f
        Extents.union(acc, _extent(f))
    end
    return _float64_xy_extent(ext)
end

function _float64_xy_extent(ext::Extents.Extent)
    xbounds = map(Float64, ext.X)
    ybounds = map(Float64, ext.Y)
    return Extents.Extent(X=xbounds, Y=ybounds)
end

