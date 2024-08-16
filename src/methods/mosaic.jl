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
function mosaic(f::Function, r1::RasterStackOrArray, rs::RasterStackOrArray...; kw...)
    mosaic(f, (r1, rs...); kw...)
end
mosaic(f::Function, regions; kw...) = _mosaic(f, first(regions), regions; kw...)
function _mosaic(f::Function, A1::AbstractRaster, regions;
    missingval=nokw,
    maskingval=nokw,
    filename=nothing,
    suffix=nothing,
    driver=nokw,
    options=nokw,
    force=false,
    kw...
)
    isnothing(missingval) && throw(ArgumentError("missingval cannot be `nothing` for `mosaic`"))
    maskingval = isnokw(maskingval) ? Rasters.missingval(first(regions)) : maskingval
    missingval = if isnokw(missingval)
        mv = Rasters.missingval(first(regions)) 
        isnokwornothing(mv) ? missing : mv
    else
        missingval
    end
    if !isnothing(filename) && (ismissing(missingval) || isnokwornothing(missingval))
        missingval = _type_missingval(eltype(A1))
    end
    T = Base.promote_type(typeof(missingval), Base.promote_eltype(regions...))
    dims = _mosaic(Tuple(map(DD.dims, regions)))
    l1 = first(regions)

    return create(filename, T, dims;
        name=name(l1),
        fill=missingval,
        missingval,
        maskingval,
        driver,
        options,
        force
    ) do C
        mosaic!(f, C, regions; missingval, kw...)
    end
end
function _mosaic(f::Function, ::AbstractRasterStack, regions;
    filename=nothing,
    suffix=keys(first(regions)),
    kw...
)
    # TODO make this write inside a single netcdf
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
mosaic!(f::Function, x::RasterStackOrArray, regions::RasterStackOrArray...; kw...) =
    mosaic!(f, x, regions; kw...)
function mosaic!(f::Function, A::AbstractRaster{T}, regions;
    missingval=Rasters.missingval(A),
    atol=_default_atol(T)
) where T
    isnokwornothing(missingval) && throw(ArgumentError("destination array must have a `missingval`"))
    _without_mapped_crs(A) do A1
        broadcast!(A1, DimSelectors(A1; atol)) do ds
            # Get all the regions that have this point
            ls = foldl(regions; init=()) do acc, l
                if DD.hasselection(l, ds)
                    v = l[ds...]
                    (acc..., l)
                else
                    acc
                end
            end
            values = foldl(ls; init=()) do acc, l
                v = l[ds...]
                if isnothing(Rasters.missingval(l))
                    (acc..., v)
                elseif ismissing(Rasters.missingval(l))
                    ismissing(v) ? acc : (acc..., v)
                else
                    v === Rasters.missingval(l) ? acc : (acc..., v)
                end
            end
            if length(values) === 0
                missingval
            else
                f(values)
            end
        end
    end
    return A
end
function mosaic!(f::Function, st::AbstractRasterStack, regions::Tuple; kw...)
    map(st, regions...) do A, r...
        mosaic!(f, A, r; kw...)
    end
end

_mosaic(alldims::Tuple{<:DimTuple,Vararg{<:DimTuple}}) = map(_mosaic, alldims...)
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

# These are pretty random default, but seem to work
_default_atol(T::Type{<:Float32}) = 100eps(T)
_default_atol(T::Type{<:Float64}) = 1000eps(T)
_default_atol(T::Type{<:Integer}) = T(1)
_default_atol(::Type) = nothing
