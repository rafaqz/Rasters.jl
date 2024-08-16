"""
    crop(x; to, touches=false, [geometrycolumn])
    crop(xs...; to)

Crop one or multiple [`AbstractRaster`](@ref) or [`AbstractRasterStack`](@ref) `x`
to match the size of the object `to`, or smallest of any dimensions that are shared.

`crop` is lazy, using a `view` into the object rather than allocating new memory.

# Keywords

- `to`: the object to crop to. This can be $OBJ_ARGUMENT
  If no `to` keyword is passed, the smallest shared area of all `xs` is used.
- `touches`: `true` or `false`. Whether to use `Touches` wraper on the object extent.

   When lines need to be included in e.g. zonal statistics, `true` should be used.
$GEOMETRYCOLUMN_KEYWORD

As `crop` is lazy, `filename` and `suffix` keywords are not used.

# Example

Crop to another raster:

```jldoctest
using Rasters, RasterDataSources, Plots
evenness = Raster(EarthEnv{HabitatHeterogeneity}, :evenness)
rnge = Raster(EarthEnv{HabitatHeterogeneity}, :range)

# Roughly cut out New Zealand from the evenness raster
nz_bounds = X(165 .. 180), Y(-50 .. -32)
nz_evenness = evenness[nz_bounds...]

# Crop range to match evenness
nz_range = crop(rnge; to=nz_evenness)
plot(nz_range)

savefig("build/nz_crop_example.png")
nothing

# output
```

![new zealand evenness cropped](nz_crop_example.png)

Crop to a polygon:

```jldoctest
using Rasters, RasterDataSources, Plots, Dates, Shapefile, Downloads

# Download a borders shapefile
shapefile_url = "https://github.com/nvkelso/natural-earth-vector/raw/master/10m_cultural/ne_10m_admin_0_countries.shp"
shapefile_name = "boundary.shp"
isfile(shapefile_name) || Downloads.download(shapefile_url, shapefile_name)
shp = Shapefile.Handle(shapefile_name).shapes[6]

evenness = Raster(EarthEnv{HabitatHeterogeneity}, :evenness)
argentina_evenness = crop(evenness; to=shp)
plot(argentina_evenness)

savefig("build/argentina_crop_example.png"); nothing

# output
```

![argentina evenness cropped](argentina_crop_example.png)

$EXPERIMENTAL
"""

function crop end
function crop(l1::RasterStackOrArray, l2::RasterStackOrArray, ls::RasterStackOrArray...; kw...)
    crop((l1, l2, ls...); kw...)
end
function crop(xs; to=nothing, kw...)
    if isnothing(to)
        to = _subsetbounds((max, min), xs)
        map(l -> _crop_to(l, to), xs)
    else
        map(l -> crop(l; to, kw...), xs)
    end
end
crop(x::RasterStackOrArray; to, geometrycolumn=nothing, kw...) = _crop_to(x, to; geometrycolumn, kw...)

# crop `A` to values of dims of `to`
function _crop_to(x, to; geometrycolumn, kw...)
    ext = _extent(to; geometrycolumn)
    if isnothing(ext)
        if isnothing(dims(to))
            throw(ArgumentError("No dims or extent available on `to` object of type $(typeof(to))"))
        else
            return _crop_to(x, dims(x, dims(to)); kw...)
        end
    else
        return _crop_to(x, ext; kw...)
    end
end
_crop_to(A, to::RasterStackOrArray; dims=DD.dims(to), kw...) = _crop_to(A, DD.dims(to, dims); kw...)
_crop_to(x, to::Dimension; kw...) = _crop_to(x, (to,); kw...)
function _crop_to(x, to::DimTuple; touches=false, kw...)
    # We can only crop to sampled dims (e.g. not categorical dims like Band)
    sampled = reduce(to; init=()) do acc, d
        lookup(d) isa AbstractSampled ? (acc..., d) : acc
    end
    return _crop_to(x, Extents.extent(sampled); touches)
end
function _crop_to(x, to::Extents.Extent; touches=false, kw...)
    # Take a view over the bounds
    _without_mapped_crs(x) do x1
        if touches
            view(x1, Touches(to))
        else
            view(x1, to)
        end
    end
end

"""
    extend(xs...; [to])
    extend(xs; [to])
    extend(x::Union{AbstractRaster,AbstractRasterStack}; to, kw...)

Extend one or multiple [`AbstractRaster`](@ref) to match the area
covered by all `xs`, or by the keyword argument `to`.

# Keywords

- `to`: the Raster or dims to extend to. If no `to` keyword is passed, the largest
    shared area of all `xs` is used.
- `touches`: `true` or `false`. Whether to use `Touches` wrapper on the object extent.
   When lines need to be included in e.g. zonal statistics, `true` shoudle be used.
$FILENAME_KEYWORD
$SUFFIX_KEYWORD

```jldoctest
using Rasters, RasterDataSources, Plots
evenness = Raster(EarthEnv{HabitatHeterogeneity}, :evenness)
rnge = Raster(EarthEnv{HabitatHeterogeneity}, :range)

# Roughly cut out South America
sa_bounds = X(-88 .. -32), Y(-57 .. 13)
sa_evenness = evenness[sa_bounds...]

# Extend range to match the whole-world raster
sa_range = extend(sa_evenness; to=rnge)
plot(sa_range)

savefig("build/extend_example.png")
nothing
# output
```

![extend](extend_example.png)

$EXPERIMENTAL
"""
function extend end
function extend(l1::RasterStackOrArray, l2::RasterStackOrArray, ls::RasterStackOrArray...; kw...)
    extend((l1, l2, ls...); kw...)
end
function extend(xs; to=nothing, kw...)
    if isnothing(to)
        to = _subsetbounds((min, max), xs)
        map(l -> _extend_to(l, to, kw...), xs)
    else
        map(l -> extend(l; to, kw...), xs)
    end
end
extend(x::RasterStackOrArray; to=dims(x), kw...) = _extend_to(x, to; kw...)

_extend_to(x::RasterStackOrArray, to::RasterStackOrArray; kw...) = _extend_to(x, dims(to); kw...)
function _extend_to(x::RasterStackOrArray, to; kw...)
    ext = _extent(to)
    isnothing(ext) && throw(ArgumentError("No dims or extent available on `to` object of type $(typeof(to))"))
    return _extend_to(x, ext; kw...)
end
_extend_to(x::RasterStackOrArray, to::Dimension; kw...) = _extend_to(x, (to,); kw...)

function _extend_to(A::AbstractRaster, to::DimTuple;
    filename=nothing,
    missingval=(isnothing(missingval(A)) ? nokw : missingval(A)),
    fill=nokw,
    touches=false,
    verbose=true,
    kw...
)
    others = otherdims(to, A)
    # Allow not specifying all dimensions
    to = (set(dims(A), map(=>, dims(A, to), to)...)..., others...)
    # Calculate the range of the old array in the extended array
    rangedims = _without_mapped_crs(A) do A
        _without_mapped_crs(to) do to
            map(dims(A, to), to) do d, t
                # Values must match exactly, so use `At`
                DD.selectindices(t, At(first(d))):DD.selectindices(t, At(last(d)))
            end
        end
    end
    others1 = otherdims(to, A)
    final_to = (set(dims(A), map(=>, dims(A, to), to)...)..., others1...)

    # If we are writing to disk swap missingval to something writeable
    if ismissing(missingval) 
        missingval = _writeable_missing(filename, eltype(A); verbose)
    end
    # If no `fill` is passed use `missingval` or zero
    if isnokw(fill)
        fill = isnokwornothing(missingval) ? zero(Missings.nonmissingtype(eltype(A))) : missingval
    end
    # Create a new extended array
    return create(filename, eltype(A), final_to;
        parent=parent(A), 
        missingval,
        name=name(A), 
        metadata=metadata(A), 
        maskingval=Rasters.missingval(A),
        verbose,
        fill,
        kw...
    ) do C
        # Input checks
        map(dims(A, to), dims(C, to)) do d1, d2
            if lookup(d1) isa Union{AbstractSampled,NoLookup}
                b1, b2 = bounds(d1), bounds(d2)
                b1[1] >= b2[1] || throw(ArgumentError("Lower bound of $(basetypeof(d1)) lookup of `$(b2[1])` are not larger than the original `$(b1[1])`"))
                b1[2] <= b2[2] || throw(ArgumentError("Upper bound of $(basetypeof(d2)) lookup of `$(b2[2])` is not larger than the original `$(b1[2])`"))
            elseif lookup(d1) isa Categorical
                map(lookup(d1)) do x
                    x in d2 || throw(ArgumentError("category $x not in new dimension"))
                end
            end
        end
        # The missingval may have changed for disk-based arrays
        if !isequal(Rasters.missingval(A), Rasters.missingval(C))
            A = replace_missing(A, Rasters.missingval(C))
        end
        # Somehow this is slow from disk?
        broadcast_dims!(identity, view(C, rangedims...), A)
    end
end
function _extend_to(st::AbstractRasterStack, to::DimTuple; suffix=keys(st), kw...)
    mapargs((A, s) -> _extend_to(A, to; suffix=s, kw...), st, suffix)
end
function _extend_to(ser::RasterSeries, to::DimTuple; kw...)
    map(x -> _extend_to(x, to; kw...), ser)
end
function _extend_to(x::RasterStackOrArray, extent::Extents.Extent{K}; kw...) where K
    shareddims = dims(x, dims(extent))
    bnds = map(val, dims(extent, shareddims))
    newdims = map(shareddims, bnds) do d, b
        l = lookup(d)
        # Use ranges for math because they have TwicePrecision magic
        # Define a range down to the lowest value,
        # but anchored at the existing value
        fl = LA.ordered_first(l); ll = LA.ordered_last(l)
        fb = first(b); lb = last(b)
        s = step(l)
        lowerrange = fb < first(bounds(l)) ? (fl:-s:fb) : (fl:-s:fl)
        upperrange = lb > last(bounds(l))  ? (ll: s:lb) : (ll: s:ll)
        if DD.order(l) isa ForwardOrdered
            newrange = last(lowerrange):s:last(upperrange)
        elseif order(d) isa ReverseOrdered
            newrange = last(upperrange):s:last(lowerrange)
        end
        newlookup = rebuild(l; data=newrange)
        rebuild(d, newlookup)
    end
    return _extend_to(x, newdims; kw...)
end

# Shared utils

# Get the largest or smallest dimensions in a tuple of AbstractRaster
function _subsetbounds(fs, layers)
    # Combine the dimensions of all layers
    dims = DD.combinedims(layers...; check=false)
    # Search through all the dimensions choosing the shortest
    alldims = map(DD.dims, layers)
    bounds = reduce(dims; init=()) do acc, d
        all(map(l -> hasdim(l, d), layers)) || return acc
        matchingdims = map(ds -> DD.dims(ds, (d,)), alldims)
        bounds = reduce(matchingdims) do a, b
            _choosebounds(fs, a, b)
        end
        return (acc..., rebuild(d, bounds))
    end
    return Extents.Extent{name(bounds)}(map(val, bounds))
end

# Choose bounds from either missing dimension
# (empty Tuple) or a comparison between two 1-Tuples
_choosebounds((f1, f2), (a,)::Tuple{<:Dimension}, ::Tuple{}) = bounds(a)
_choosebounds((f1, f2), (a,)::Tuple{<:Dimension}, b::Tuple{<:Dimension}) = _choosebounds((f1, f2), bounds(a), b)
_choosebounds((f1, f2), ::Tuple{}, ::Tuple{}) = ()
_choosebounds((f1, f2), ::Tuple{}, (b,)::DimTuple) = bounds(b)
_choosebounds((f1, f2), a::Tuple, ::Tuple{}) = a
function _choosebounds((f1, f2), (a1, a2)::Tuple{<:Any,<:Any}, (b,)::Tuple{<:Dimension})
    b1, b2 = bounds(b)
    return f1(a1, b1), f2(a2, b2)
end
