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
using Rasters, RasterDataSources, ArchGDAL, Plots
st = RasterStack(WorldClim{Climate}; month=1);

africa = st[X(-20.0 .. 60.0), Y(-40.0 .. 35.0)]
a = plot(africa)

aus = st[X(100.0 .. 160.0), Y(-50.0 .. -10.0)]
b = plot(aus)

# Combine with mosaic
mos = mosaic(first, aus, africa)
c = plot(mos)

savefig(a, "build/mosaic_example_africa.png")
savefig(b, "build/mosaic_example_aus.png")
savefig(c, "build/mosaic_example_combined.png")
nothing
# output

```

### Individual continents

![arica](mosaic_example_africa.png)

![aus](mosaic_example_aus.png)

### Mosaic of continents

![mosaic](mosaic_example_combined.png)

$EXPERIMENTAL
"""
function mosaic(last, r1::RasterStackOrArray, rs::RasterStackOrArray...; kw...)
    mosaic(f, (r1, rs...); kw...)
end
function mosaic(f::Function, r1::RasterStackOrArray, rs::RasterStackOrArray...; kw...)
    mosaic(f, (r1, rs...); kw...)
end
mosaic(f::Function, regions; kw...) = _mosaic(f, first(regions), regions; kw...)
function _mosaic(f::Function, ::AbstractRaster, regions;
    missingval=missingval(first(regions)), filename=nothing, suffix=nothing, kw...
)
    missingval = missingval isa Nothing ? missing : missingval
    zeros = map(r -> zero(Missings.nonmissingtype(eltype(r))), regions)
    T = Base.promote_type(typeof(missingval), typeof(f(zeros)))
    dims = _mosaic(Tuple(map(DD.dims, regions)))
    l1 = first(regions)
    A = create(filename, T, dims; name=name(l1), missingval, metadata=metadata(l1))
    open(A; write=true) do a
        mosaic!(f, a, regions; missingval, kw...)
    end
    return A
end
function _mosaic(f::Function, ::AbstractRasterStack, regions;
    filename=nothing, suffix=keys(first(regions)), kw...
)
    layers = map(suffix, map(values, regions)...) do s, A...
        mosaic(f, A...; filename, suffix=s, kw...)
    end
    return DD.rebuild_from_arrays(first(regions), Tuple(layers))
end

"""
    mosaic!(f, x, regions...; missingval, atol)
    mosaic!(f, x, regions::Tuple; missingval, atol)

Combine `regions` in `x` using the function `f`.

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

Cut out Australia and Africa stacks, then combined them
into a single stack.

```jldoctest
using Rasters, RasterDataSources, ArchGDAL, Statistics, Plots
st = read(RasterStack(WorldClim{Climate}; month=1))
aus = st[X=100.0 .. 160.0, Y=-50.0 .. -10.0]
africa = st[X=-20.0 .. 60.0, Y=-40.0 .. 35.0]
mosaic!(first, st, aus, africa)
plot(st)
savefig("build/mosaic_bang_example.png")
nothing
# output

```

![mosaic](mosaic_bang_example.png)

$EXPERIMENTAL
"""
mosaic!(dest::RasterStackOrArray, regions::RasterStackOrArray...; kw...) =
    mosaic!(last, dest, regions; kw...)
mosaic!(f::Function, dest::RasterStackOrArray, regions::RasterStackOrArray...; kw...) =
    mosaic!(f, dest, regions; kw...)
mosaic!(f, dest::RasterStackOrArray, regions::RasterStackOrArray...; kw...) =
    mosaic!(f, dest, regions; kw...)
function mosaic!(f::Function, dest::AbstractRaster, regions; kw...)
    any_intersect = any(enumerate(regions)) do (i, r)
        any(i:length(regions)) do i
            Extents.intersects(extent(r), extent(regions[i]))
        end
    end
    if any_intersect
        _mosaic!(f, dest, regions; kw...)
    else
        _no_ovelap_mosaic!(f, dest, regions; kw...)
    end
end

_mosaic!(f::typeof(first), dest::AbstractRaster, regions; kw...) =
    _mosaic!(last, dest, reverse(regions); kw...)
function _mosaic!(f::typeof(last), dest::AbstractRaster, regions;
    missingval=missingval(A), selectors=Near()
)
    # For `last` we write each region in sequence
    map(regions) do r
        # View the dest with the extent of the region
        v = view(dest, extent(r))
        # Then use nearest neighbor selector indexing to update it
        v .= parent(view(r, DimSelectors(v; selectors)))
    end

    return dest
end
function _mosaic!(f::typeof(sum), dest::AbstractRaster, regions;
    missingval=missingval(A), selectors=Near()
)
    # For `sum` we first zero the mosaic regions of dest
    map(regions) do r
        # View the dest of the extent with zeros
        view(dest, extent(r)) .+= zero(eltype(dest))
    end
    # Then sum into them
    map(regions) do r
        # Use nearest neighbor selector indexing to update dest
        v = view(dest, extent(r))
        v .+= parent(view(r, DimSelectors(v; selectors)))
    end
    return dest
end
function _mosaic!(f::typeof(mean), dest::AbstractRaster, regions; kw...)
    # First sum
    _mosaic!(sum, dest, regions; kw...)

    # Then count
    counts = similar(dest, UInt8)
    fill!(counts, 0x00)
    map(regions) do r
        view(dest, extent(r)) .+= 0x01
    end
    # Then divide the sums by the counts
    dest ./= counts

    return dest
end
function _mosaic!(f::Function, st::AbstractRasterStack, regions::Tuple; 
    dims=(Val{XDim}(), Val{YDim}()),
    pixel_done=falses(size(st, dims)),
    region_done=falses(length(regions)),
    region_intersects=falses(length(regions)),
    kw...
)
    map(st, regions...) do A, r...
        pixel_done .= false
        region_done .= false
        mosaic!(f, A, r; dims, pixel_done, region_done, region_intersects, kw...)
    end
end
function _mosaic!(f::Function, A::AbstractRaster, regions;
    missingval=missingval(A),
    selectors=Near(),
    dims=(Val{XDim}(), Val{YDim}()),
    pixel_done=falses(size(A, dims)),
    region_done=falses(length(regions)),
    region_intersects=falses(length(regions)),
)
    region_vect = collect(regions)
    for (i1, r1) in enumerate(regions)
        region_intersects .= false
        region_done[i1] && continue
        for (i2, r2) in  enumerate(regions)
            region_done[i2] && continue
            region_intersects[i2] = intersects(extent(r1), extent(r2))
            if contains(extent(r1), extent(r2))
                region_done[i2] = true
            end
        end
        intersecting = region_vec[region_intersects]
        _mosaic_intersecting!(f, dest, pixel_done, intersecting)
    end
end

function _mosaic_intersecting!(f, dest, intersecting)
    n = length(intersecting)
    ext = mapreduce(extent, Extents.intersection, view(intersecting, j:n))
end

# Nothing overlaps, just write
function _no_ovelap_mosaic!(f, dest, regions; kw...)
    g(x) = f((x,))
    map(regions) do r
        # Use nearest neighbor selector indexing to update dest
        v = view(dest, extent(r))
        v += g.(parent(view(r, DimSelectors(v; selectors))))
    end
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
