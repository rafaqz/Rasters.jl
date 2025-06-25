const RasterVecOrTuple = Union{Tuple{Vararg{AbstractRaster}},AbstractArray{<:AbstractRaster}}
const RasterStackVecOrTuple = Union{Tuple{Vararg{AbstractRasterStack}},AbstractArray{<:AbstractRasterStack}}

const MOSAIC_ARGUMENTS = """
- `f`: A reducing function for values where `regions` overlap. 
    Note that common base functions 
    (`mean`, `sum`, `prod`, `first`, `last`, `minimum`, `maximum`,  `length`)
    are optimised and will work on many memory or disk based files,
    but user-defined functions may fail at larger scales unless `op` is passes as a keyword.
- `regions`: Iterable of `Raster` or `RasterStack`. Using an `AbstractArray` is 
    usually better than a `Tuple` or splat when there are many regions.
"""

const MOSAIC_KEYWORDS = """
- `missingval`: Fills empty areas, and defaults to the
    `missingval` of the first region.
- `op`: an operator for the reduction, e.g. `add_sum` for `sum`. 
    For common methods like `sum` these are known and detected for you, 
    but you can provide it manually for other functions, so they continue
    to work at large scales.
- `atol`: Absolute tolerance for comparison between index values.
    This is often required due to minor differences in range values
    due to floating point error. It is not applied to non-float dimensions.
    A tuple of tolerances may be passed, matching the dimension order.
$PROGRESS_KEYWORD
- `read`: `read` lazy raster regions before writing them. This may help
    if there are chunk alignment issues between the source and dest rasters.
    `false` by default.
"""

"""
    mosaic(f, regions...; kw...)
    mosaic(f, regions; kw...)

Combine `regions` of `Raster` or `RasterStack` into a single object, 
using the function `f` to combine overlapping areas.

# Arguments

$MOSAIC_ARGUMENTS

# Keywords

$MOSAIC_KEYWORDS
$FILENAME_KEYWORD
$SUFFIX_KEYWORD
$FORCE_KEYWORD

If your mosaic has has apparent line errors, increase the `atol` value.

# Example

Here we cut out Australia and Africa from a stack, and join them with `mosaic`.

```jldoctest
using Rasters, RasterDataSources, NaturalEarth, DataFrames, Dates, Plots
import ArchGDAL

countries = naturalearth("admin_0_countries", 110) |> DataFrame
climate = RasterStack(WorldClim{Climate}, (:tmin, :tmax, :prec, :wind); month=July)
country_climates = map(("Norway", "Denmark", "Sweden")) do name
    country = subset(countries, :NAME => ByRow(==(name)))
    trim(mask(climate; with=country); pad=10)
end
scandinavia_climate = trim(mosaic(first, country_climates; progress=false))
plot(scandinavia_climate)
savefig("build/mosaic_example_combined.png"); nothing

# output

```

### Mosaic of countries

![mosaic](mosaic_example_combined.png)

$EXPERIMENTAL
"""
mosaic(f::Function, r1::RasterStackOrArray, rs::RasterStackOrArray...; kw...) =
    _mosaic(f, r1, [r1, rs...]; kw...)
mosaic(f::Function, regions; kw...) = _mosaic(f, first(regions), regions; kw...)
_mosaic(f::Function, R1::RasterStackOrArray, regions; kw...) = 
    _mosaic(f, R1, collect(regions); kw...)
function _mosaic(f::Function, R1::RasterStackOrArray, regions::AbstractArray;
    to=nothing,
    filename=nothing,
    atol=nothing,
    missingval=nokw,
    op=nokw,
    progress=true,
    read=false,
    kw...
)
    dims = if isnothing(to)
        _mosaic(map(DD.dims, regions))
    else
        ds = DD.dims(to)
        isnothing(ds) && throw(ArgumentError("`to` object does not return Dimensions from `dims(to)`. Pass a tuple of Dimension, a Raster, or RasterStack"))
        ds
    end

    missingval = if isnothing(missingval) 
        throw(ArgumentError("missingval cannot be `nothing` for `mosaic`"))
    elseif isnokw(missingval)
        mv = Rasters.missingval(first(regions)) 
        isnokwornothing(mv) ? missing : mv
    else
        missingval
    end
    missingval_pair = if missingval isa Pair
        missingval
    elseif !isnothing(filename) && (ismissing(missingval) || isnokw(missingval))
        _type_missingval(eltype(R1)) => missing
    else
        missingval => missingval
    end
    T = _mosaic_eltype(f, R1, regions)

    return create(filename, T, dims;
        name=name(R1),
        fill=missingval_pair[1],
        missingval=missingval_pair,
        kw...
    ) do C
        mosaic!(f, C, regions; op, atol, progress, read)
    end
end

function _mosaic_eltype(f, R1::AbstractRaster, regions)
    E = reduce(regions; init=eltype(R1)) do T, r
        promote_type(T, eltype(r))
    end
    V = Vector{E}
    return Base.promote_op(f, V)
end
function _mosaic_eltype(f, R1::AbstractRasterStack{K}, regions) where K
    map(K) do k
        E = reduce(regions; init=eltype(R1[k])) do T, r
            promote_type(T, eltype(r[k]))
        end
        V = Vector{E}
        Base.promote_op(f, V)
    end |> NamedTuple{K}
end

"""
    mosaic!(f, dest::Union{Raster,RasterStack}, regions...; kw...)
    mosaic!(f, dest::Union{Raster,RasterStack}, regions; kw...)

Combine `regions` of `Raster` or `RasterStack` into `dest`
using the function `f` to combine overlapping areas.

# Arguments

$MOSAIC_ARGUMENTS

- `dest`: A `Raster` or `RasterStack`. May be a an opened disk-based `Raster`,
    the result will be written to disk.
    With the current algorithm, the read speed is slow.

# Keywords

$MOSAIC_KEYWORDS

$EXPERIMENTAL
"""
mosaic!(f::Function, dest::RasterStackOrArray, regions::RasterStackOrArray...; kw...) =
    mosaic!(f, dest, regions; kw...)
function mosaic!(
    f::Function, 
    dest::RasterStackOrArray,
    regions::Union{Tuple,AbstractArray}; 
    op=nokw,
    kw...
)
    op = isnokw(op) ? _reduce_op(f, missingval(dest)) : op
    # Centering avoids pixel edge floating point error
    dest_centered = _prepare_for_burning(dest; order=nothing)
    regions_centered = map(r -> _prepare_for_burning(r; order=nothing), regions)
    _mosaic!(f, op, dest_centered, regions_centered; kw...)
    return dest
end

function _mosaic!(
    f::typeof(mean), op::Nothing, dest::Raster, regions::RasterVecOrTuple;
    kw...
)
    if length(regions) <= typemax(UInt8)
        _mosaic_mean!(dest, UInt8, regions; kw...)
    elseif length(regions) <= typemax(UInt16)
        _mosaic_mean!(dest, UInt16, regions; kw...)
    else
        _mosaic_mean!(dest, UInt32, regions; kw...)
    end
end
function _mosaic!(
    f::typeof(length), op::Nothing, dest::AbstractRaster, regions::RasterVecOrTuple;
    read=false, kw...
)
    for region in regions
        if read
            _count_region!(dest, region; kw...)
        else
            _count_region!(dest, Base.read(region); kw...)
        end
    end
    return dest
end
# Where there is a known reduction operator we can apply each region as a whole
function _mosaic!(
    f::Function, op, dest::AbstractRaster, regions::RasterVecOrTuple; 
    gc::Union{Integer,Nothing}=nothing,
    read=false,
    progress=true,
    _progressmeter=_mosaic_progress(f, progress, length(regions)),
    kw...
)
    for (i, region) in enumerate(regions)
        !isnothing(_progressmeter) && ProgressMeter.next!(_progressmeter)
        # Run garbage collector every `gc` regions
        !isnothing(gc) && rem(i, gc) == 0 && GC.gc()
        if read
            _mosaic_region!(op, dest, Base.read(region); kw...)
        else
            open(region) do R
                _mosaic_region!(op, dest, R; kw...)
            end
        end
    end
    return dest
end
# Generic unknown functions
function _mosaic!(
    f::Function, op::Nothing, A::AbstractRaster{T}, regions::RasterVecOrTuple;
    missingval=missingval(A), 
    atol=nothing,
    progress=false,
    read=false,
) where T
    # TODO: rethink this algorithm using overlapping subregions
    # Its super slow to do this for every pixel
    isnokwornothing(missingval) && throw(ArgumentError("destination array must have a `missingval`"))
    R = promote_type(map(Missings.nonmissingtype ∘ eltype, regions)...)
    buffer = Vector{R}(undef, length(regions))
    _without_mapped_crs(A) do A1
        broadcast!(A1, DimSelectors(A1; atol)) do ds
            # Get all the regions that have this point
            i = 0
            for r in regions 
                if DD.hasselection(r, ds)
                    x = r[ds...]
                    if x !== Rasters.missingval(r)
                        i += 1
                        buffer[i] = x
                    end
                end
            end
            if i === 0
                missingval
            else
                f(view(buffer, 1:i))
            end
        end
    end
    return A
end
function _mosaic!(
    f::Function, 
    op, 
    st::AbstractRasterStack{K}, 
    regions::RasterStackVecOrTuple; 
    progress=true,
    kw...
) where K
    _progressmeter = _mosaic_progress(f, progress, length(regions), length(K))
    map(values(st), K) do A, k
        layer_regions = map(r -> r[k], regions)
        _mosaic!(f, op, A, layer_regions; kw..., _progressmeter)
    end
    return st
end

function _mosaic_mean!(dest, ::Type{T}, regions; 
    read=false, kw...
) where T
    # Note: sum and count are separate broadcasts because 
    # most disk formats don't support writing a tuple

    # Define a Raster to count into
    counts = create(nothing, T, dest; missingval=zero(T))
    counts .= zero(T)
    for region in regions
        # Add region to dest
        if read 
            region1 = Base.read(region)
            _mosaic_region!(Base.add_sum, dest, region1; kw...)
            _count_region!(counts, region1; kw...)
        else
            _mosaic_region!(Base.add_sum, dest, region; kw...)
            _count_region!(counts, region; kw...)
        end
    end
    # Divide dest by counts
    # Avoid divide by zero for missing values
    dest .= ((d, c) -> d === missingval(dest) ? missingval(dest) : d / c).(dest, counts)
    return dest
end

function _mosaic_progress(f, progress, n, l=nothing)
    if progress
        if isnothing(l)
            _progress(n; desc="Mosaicing $n regions with $f...")
        else
            _progress(n * l; desc="Mosaicing $l layers of $n regions with $f...")
        end
    else
        nothing
    end
end

function _mosaic_region!(op, dest, region; atol=nothing, kw...)
    function skip_or_op(a, b) 
        if b === missingval(region)
            a
        elseif a === missingval(dest) 
            b
        else
            op(a, b)
        end
    end
    ext = _maybe_pad_floats(extent(region), sampling(dest))
    selectors = map(sampling(dest)) do sa
        ispoints(sa) ?  At(; atol) : Contains()
    end
    ds = DimSelectors(view(dest, ext); selectors)
    # `parent` needed to skip broadcast checks
    dest[ext] .= skip_or_op.(parent(view(dest, ext)), parent(view(region, ds)))
    return dest
end

function _count_region!(count::AbstractRaster{T}, region::AbstractRaster; kw...) where T
    function skip_or_count(a, b)
        if b === missingval(region)
            a
        elseif a === missingval(count) 
            oneunit(Missings.nonmissingtype(T))
        else
            a + oneunit(a)
        end
    end
    ext = extent(region)
    ds = DimSelectors(view(count, ext))
    view(count, ext) .= skip_or_count.(view(count, ext), view(region, ds))
    return count
end


function _mosaic(dimtuplevec::AbstractArray{<:DimTuple})
    map(first(dimtuplevec)) do d
        dimvec = map(ds -> dims(ds, d), dimtuplevec)
        _mosaic(dimvec)
    end
end
function _mosaic(dimsvec::AbstractArray{<:Dimension})
    d1 = first(dimsvec)
    map(dimsvec) do d
        DD.comparedims(d1, d; val=false, length=false, valtype=true)
    end
    return rebuild(d1, _mosaic(map(lookup, dimsvec)))
end
_mosaic(lookups::AbstractArray{<:Lookup}) = _mosaic(first(lookups), lookups)
function _mosaic(lookup::Categorical, lookups::AbstractArray{<:Lookup})
    newindex = union(lookups...)
    if order isa ForwardOrdered
        newindex = sort(newindex; order=LA.ordering(order(lookup)))
    end
    return rebuild(lookup; data=newindex)
end
function _mosaic(lookup::AbstractSampled, lookups::AbstractArray{<:Lookup})
    order(lookup) isa Unordered && throw(ArgumentError("Cant mozaic an Unordered lookup"))
    return _mosaic(span(lookup), lookup, lookups)
end
function _mosaic(span::Regular, lookup::AbstractSampled, lookups::AbstractArray{<:Lookup})
    mi = minimum(first, lookups)
    ma = maximum(last, lookups)
    newindex = if order(lookup) isa ForwardOrdered
        if mi isa AbstractFloat
            # Handle slight range errors to make sure
            # we don't drop one step of the range
            mi:step(span):ma + 2eps(ma)
        else
            mi:step(span):ma
        end
    else
        if mi isa AbstractFloat
            ma:step(span):mi - 2eps(mi)
        else
            ma:step(span):mi
        end
    end
    return rebuild(lookup; data=newindex)
end
function _mosaic(::Irregular, lookup::AbstractSampled, lookups::AbstractArray{<:Lookup})
    newindex = sort(union(map(parent, lookups)...); order=LA.ordering(order(lookup)))
    return rebuild(lookup; data=newindex)
end
function _mosaic(span::Explicit, lookup::AbstractSampled, lookups::AbstractArray{<:Lookup})
    # TODO make this less fragile to floating point inaccuracy
    newindex = sort(union(map(parent, lookups)...); order=LA.ordering(order(lookup)))
    bounds = map(val ∘ DD.span, lookups)
    lower = map(b -> view(b, 1, :), bounds)
    upper = map(b -> view(b, 2, :), bounds)
    newlower = sort(union(lower...); order=LA.ordering(order(lookup)))
    newupper = sort(union(upper...); order=LA.ordering(order(lookup)))
    newbounds = vcat(permutedims(newlower), permutedims(newupper))
    return rebuild(lookup; data=newindex, span=Explicit(newbounds))
end

# Pad floats for intervals so that small floating point 
# error doesn't exclude values in nearly matching lookups
function _maybe_pad_floats(ext::Extent{K}, sampling::Tuple) where K
    map(values(Extents.bounds(ext)), sampling) do b, sa
        if isintervals(sa) && eltype(first(b)) <: AbstractFloat
            b[1] - 10eps(b[1]), b[2] + 10eps(b[2])
        else
            b
        end
    end |> Extent{K}
end