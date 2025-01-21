const RasterVecOrTuple = Union{Tuple{Vararg{AbstractRaster}},AbstractArray{<:AbstractRaster}}
const RasterStackVecOrTuple = Union{Tuple{Vararg{AbstractRasterStack}},AbstractArray{<:AbstractRasterStack}}

"""
    mosaic(f, regions...; missingval, atol)
    mosaic(f, regions; missingval, atol)

Combine `regions` into a single raster.

# Arguments

- `f`: A reducing function (`mean`, `sum`, `first`, `last` etc.)
    for values where `regions` overlap.
- `regions`: Iterable or splatted `Raster` or `RasterStack`.

# Keywords

- `missingval`: Fills empty areas, and defualts to the
    `missingval` of the first region.
- `atol`: Absolute tolerance for comparison between index values.
    This is often required due to minor differences in range values
    due to floating point error. It is not applied to non-float dimensions.
    A tuple of tolerances may be passed, matching the dimension order.
$FILENAME_KEYWORD
$SUFFIX_KEYWORD

If your mosaic has has apparent line errors, increase the `atol` value.

# Example

Here we cut out Australia and Africa from a stack, and join them with `mosaic`.

```jldoctest
using Rasters, RasterDataSources, NaturalEarth, DataFrames, Dates
countries = naturalearth("admin_0_countries", 110) |> DataFrame
climate = RasterStack(WorldClim{Climate}, (:tmin, :tmax, :prec, :wind); month=July)
country_climates = map(("Norway", "Denmark", "Sweden")) do name
    country = subset(countries, :NAME => ByRow(==("Norway")))
    trim(mask(climate; with=country); pad=10)
end
scandinavia_climate = trim(mosaic(first, country_climates))
plot(scandinavia_climate)

savefig("build/mosaic_example_combined.png")
# output

```

### Mosaic of countries

![mosaic](mosaic_example_combined.png)

$EXPERIMENTAL
"""
mosaic(f::Function, r1::RasterStackOrArray, rs::RasterStackOrArray...; kw...) =
    mosaic(f, (r1, rs...); kw...)
mosaic(f::Function, regions; kw...) = _mosaic(f, first(regions), regions; kw...)
function _mosaic(f::Function, A1::AbstractRaster, regions;
    missingval=nokw,
    filename=nothing,
    suffix=nothing,
    driver=nokw,
    options=nokw,
    force=false,
    kw...
)
    V = Vector{promote_type(map(Missings.nonmissingtype ∘ eltype, regions)...)}
    T = Base.promote_op(f, V)
    dims = _mosaic(Tuple(map(DD.dims, regions)))
    l1 = first(regions)

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
        _type_missingval(eltype(A1)) => missing
    else
        missingval => missingval
    end

    return create(filename, T, dims;
        name=name(l1),
        fill=missingval_pair[1],
        missingval=missingval_pair,
        driver,
        options,
        force
    ) do C
        mosaic!(f, C, regions; missingval, kw...)
    end
end
function _mosaic(f::Function, r1::AbstractRasterStack, regions;
    filename=nothing,
    suffix=keys(first(regions)),
    kw...
)
    # TODO make this write inside a single netcdf
    layers = map(suffix, map(values, regions)...) do s, A...
        mosaic(f, A...; filename, suffix=s, kw...)
    end
    return DD.rebuild_from_arrays(r1, layers)
end

"""
    mosaic!(f, x, regions...; missingval, atol)
    mosaic!(f, x, regions::Tuple; missingval, atol)

Combine `regions` in `Raster` or `RasterStack` `x` using the function `f`.

# Arguments

- `f` a function (e.g. `mean`, `sum`, `first` or `last`) that is applied to
    values where `regions` overlap.
- `x`: A `Raster` or `RasterStack`. May be a an opened disk-based `Raster`,
    the result will be written to disk.
    With the current algorithm, the read speed is slow.
- `regions`: source objects to be joined. These should be memory-backed
    (use `read` first), or may experience poor performance. If all objects have
    the same extent, `mosaic` is simply a merge.

# Keywords

- `missingval`: Fills empty areas, and defualts to the `missingval/
    of the first layer.
- `atol`: Absolute tolerance for comparison between index values.
    This is often required due to minor differences in range values
    due to floating point error. It is not applied to non-float dimensions.
    A tuple of tolerances may be passed, matching the dimension order.

# Example

Cut out scandinavian countries and plot:

```jldoctest
using Rasters, RasterDataSources, NaturalEarth, DataFrames, Dates
# Get climate data form worldclim
climate = RasterStack(WorldClim{Climate}, (:tmin, :tmax, :prec, :wind); month=July)
# And country borders from natural earth
countries = naturalearth("admin_0_countries", 110) |> DataFrame
# Cut out each country
country_climates = map(("Norway", "Denmark", "Sweden")) do name
    country = subset(countries, :NAME => ByRow(==(name)))
    trim(mask(climate; with=country); pad=10)
end
# Mosaic together to a single raster
scandinavia_climate = mosaic(first, country_climates)
# And plot
plot(scandinavia_climate)

savefig("build/mosaic_bang_example.png")
# output

```

![mosaic](mosaic_bang_example.png)

$EXPERIMENTAL
"""
mosaic!(f::Function, dest::RasterStackOrArray, regions::RasterStackOrArray...; kw...) =
    mosaic!(f, dest, regions; kw...)
function mosaic!(
    f::Function, 
    dest::RasterStackOrArray,
    regions::Union{Tuple,AbstractArray}; 
    op=_reduce_op(f, missingval(dest)), 
    kw...
)
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
    kw...
)
    for region in regions
        _count_region!(dest, region; kw...)
    end
    return dest
end
# Where there is a known reduction operator we can apply each region as a whole
function _mosaic!(
    f::Function, op, dest::AbstractRaster, regions::RasterVecOrTuple; 
    kw...
)
    for region in regions
        _mosaic_region!(op, dest, region; kw...)
    end
    return dest
end
# Generic unknown functions
function _mosaic!(
    f::Function, op::Nothing, A::AbstractRaster{T}, regions::RasterVecOrTuple;
    missingval=missingval(A), atol=nothing
) where T
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
    st::AbstractRasterStack, 
    regions::RasterStackVecOrTuple; 
    kw...
)
    map(values(st), map(values, regions)...) do A, r...
        _mosaic!(f, op, A, r; kw...)
    end
    return st
end

function _mosaic_mean!(dest, ::Type{T}, regions; kw...) where T
    # Note: sum and count are separate broadcasts because 
    # most disk formats don't support writing a tuple

    # Define a Raster to count into
    counts = create(nothing, T, dest; missingval=zero(T))
    counts .= zero(T)
    for region in regions
        # Add region to dest
        _mosaic_region!(Base.add_sum, dest, region; kw...)
        # Count region
        _count_region!(counts, region; kw...)
    end
    # Divide dest by counts
    # Avoid divide by zero for missing values
    dest .= ((d, c) -> d === missingval(dest) ? missingval(dest) : d / c).(dest, counts)
    return dest
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

_mosaic(alldims::Tuple{<:DimTuple,Vararg{DimTuple}}) = map(_mosaic, alldims...)
function _mosaic(dims::Dimension...)
    map(dims) do d
        DD.comparedims(first(dims), d; val=false, length=false, valtype=true)
    end
    return rebuild(first(dims), _mosaic(lookup(dims)))
end
_mosaic(lookups::LookupTuple) = _mosaic(first(lookups), lookups)
function _mosaic(lookup::Categorical, lookups::LookupTuple)
    newindex = union(lookups...)
    if order isa ForwardOrdered
        newindex = sort(newindex; order=LA.ordering(order(lookup)))
    end
    return rebuild(lookup; data=newindex)
end
function _mosaic(lookup::AbstractSampled, lookups::LookupTuple)
    order(lookup) isa Unordered && throw(ArgumentError("Cant mozaic an Unordered lookup"))
    return _mosaic(span(lookup), lookup, lookups)
end
function _mosaic(span::Regular, lookup::AbstractSampled, lookups::LookupTuple)
    newindex = if order(lookup) isa ForwardOrdered
        mi = minimum(map(first, lookups))
        ma = maximum(map(last, lookups))
        if mi isa AbstractFloat
            # Handle slight range errors to make sure
            # we don't drop one step of the range
            mi:step(span):ma + 2eps(ma)
        else
            mi:step(span):ma
        end
    else
        mi = minimum(map(last, lookups))
        ma = maximum(map(first, lookups))
        if mi isa AbstractFloat
            ma:step(span):mi - 2eps(mi)
        else
            ma:step(span):mi
        end
    end
    return rebuild(lookup; data=newindex)
end
function _mosaic(::Irregular, lookup::AbstractSampled, lookups::LookupTuple)
    newindex = sort(union(map(parent, lookups)...); order=LA.ordering(order(lookup)))
    return rebuild(lookup; data=newindex)
end
function _mosaic(span::Explicit, lookup::AbstractSampled, lookups::LookupTuple)
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
# error doesn't exclude values in nealy matching lookups
function _maybe_pad_floats(ext::Extent{K}, sampling::Tuple) where K
    map(values(Extents.bounds(ext)), sampling) do b, sa
        if isintervals(sa) && eltype(first(b)) <: AbstractFloat
            b[1] - 10eps(b[1]), b[2] + 10eps(b[2])
        else
            b
        end
    end |> Extent{K}
end