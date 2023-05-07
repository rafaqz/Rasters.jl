_take_last(a, b) = b

_reduce_op(::typeof(sum)) = Base.add_sum
_reduce_op(::typeof(prod)) = Base.mul_prod
_reduce_op(::typeof(minimum)) = min
_reduce_op(::typeof(maximum)) = max
_reduce_op(x) = nothing

_reduce_init(reduce, st::AbstractRasterStack) = map(A -> _reduce_init(reduce, A), st)
_reduce_init(reduce, ::AbstractRaster{T}) where T = _reduce_init(reduce, T)
_reduce_init(reduce, nt::NamedTuple) = map(x -> _reduce_init(reduce, x), nt)
_reduce_init(f, x) = _reduce_init(f, typeof(x))

_reduce_init(::Nothing, x::Type{T}) where T = zero(T)
_reduce_init(f::Function, ::Type{T}) where T = zero(f((zero(nonmissingtype(T)), zero(nonmissingtype(T)))))
_reduce_init(::typeof(sum), ::Type{T}) where T = zero(nonmissingtype(T))
_reduce_init(::typeof(prod), ::Type{T}) where T = oneunit(nonmissingtype(T))
_reduce_init(::typeof(minimum), ::Type{T}) where T = typemax(nonmissingtype(T))
_reduce_init(::typeof(maximum), ::Type{T}) where T = typemin(nonmissingtype(T))


# Handles all data validations needed to
# run before rasterizing
struct Rasterizer{T,G,F,FV,FI,I,R,O,M,MD,A,ET}
    to::T
    geom::G
    fill::F
    fillval::FV
    fillitr::FI
    reduce::R
    op::O
    init::I
    name::Symbol
    missingval::M
    metadata::MD
    suffix::String
    filename::Union{String,Nothing}
    shape::Symbol
    boundary::Symbol
    verbose::Bool
    progress::Bool
    allocs::A
    lock::SectorLocks
    eltype::ET
    crs::Any
    mappedcrs::Any
end
function Rasterizer(to, geom, fill, fillitr;
    name=Symbol(""),
    reduce=nothing,
    op=_reduce_op(reduce),
    missingval=nothing,
    metadata=Metadata(Dict()),
    shape=nothing,
    eltype=nothing,
    suffix="",
    verbose=true,
    progress=true,
    filename=nothing,
    boundary=:center,
    res=nothing, # We shouldn't need this but coverage does
    crs=nothing,
    mappedcrs=nothing,
    init=nothing,
)
    shape = if isnothing(shape)
        if GI.isgeometry(geom)
            _geom_shape(geom)
        else
            _geom_shape(first(geom))
        end
    else
        shape
    end

    filleltype = if fillitr isa NamedTuple 
        if all(map(x -> x isa Number, fillitr))
            map(typeof, fillitr)
        elseif all(map(x -> Base.IteratorEltype(x) isa Base.HasEltype, fillitr))
            map(Base.eltype, fillitr)
        else
            map(typeof ∘ first, fillitr) # This is not really correct. 
        end
    elseif fill isa Function
        isnothing(init) ? typeof(missingval) : typeof(fill(init))
    elseif Base.IteratorEltype(fillitr) isa Base.HasEltype
        Base.eltype(fillitr)
    else
        typeof(first(fillitr)) # This is not really correct
    end

    stable_reductions = (first, last, sum, prod, maximum, minimum)
    init = isnothing(init) ? _reduce_init(reduce, filleltype) : init
    if shape == :points &&
        !GI.isgeometry(geom) && 
        !GI.trait(first(geom)) isa PointTrait &&
        !(reduce in stable_reductions)
        @warn "currently `:points` rasterization of multiple non-`PointTrait` geometries may be innaccurate for `reduce` methods besides $stable_reductions. Make a Rasters.jl github issue if you need this to work"
    end
    name = Symbol(_filter_name(name, fill))
    eltype, missingval = get_eltype_missingval(eltype, missingval, fillitr, init, filename, op, reduce)
    lock = SectorLocks()
    allocs = shape == :points ? nothing : _burning_allocs(to)
    # shape === :point && @warn "Rasterizing geoms as a :point, not $shape"
    to = _as_intervals(to) # Only makes sense to rasterize to intervals

    return Rasterizer(to, geom, fill, fill, fillitr, reduce, op, init, name, missingval, metadata,
                      suffix, filename, shape, boundary, verbose, progress, allocs, lock, eltype, crs, mappedcrs)
end

function get_eltype_missingval(eltype, missingval, fillitr, init::NamedTuple, filename, op, reduce)
    eltype = eltype isa NamedTuple ? eltype : map(_ -> eltype, init)
    missingval = missingval isa NamedTuple ? missingval : map(_ -> missingval, init)
    em = map(eltype, missingval, fillitr, init) do et, mv, fi, i
        get_eltype_missingval(et, mv, fi, i, filename, op, reduce)
    end
    eltype = map(first, em)
    missingval = map(last, em)
    return eltype, missingval
end
function get_eltype_missingval(known_eltype, missingval, fillitr, init, filename, op, reduce)
    eltype = if isnothing(known_eltype)
        if fillitr isa Function 
            typeof(fillitr(init))
        elseif op isa Function
            typeof(op(init, init))
        elseif reduce isa Function
            typeof(reduce((init, init)))
        else
            typeof(init)
        end
    else
        known_eltype
    end

    missingval = isnothing(missingval) ? _writeable_missing(filename, eltype) : missingval
    # eltype was not the actually array eltype, so promote it with the missingval
    eltype = isnothing(known_eltype) ? promote_type(typeof(missingval), eltype) : eltype
    return eltype, missingval
end

function Rasterizer(to::AbstractRaster, data; missingval=missingval(to), name=name(to), kw...)
    Rasterizer(dims(to), data; kw..., missingval, name)
end
function Rasterizer(to::AbstractRasterStack, data; fill, name=keys(to), kw...)
    Rasterizer(dims(to), data; fill, name, kw...)
end
function Rasterizer(to::Nothing, data; fill, kw...)
    to = _extent(data)
    Rasterizer(to, data; kw..., fill)
end
function Rasterizer(to, data; fill, kw...)
    Rasterizer(_extent(to), data; kw..., fill)
end
function Rasterizer(to::Extents.Extent{K}, data;
    fill,
    res::Union{Nothing,Real,NTuple{<:Any,<:Real}}=nothing,
    size::Union{Nothing,Int,NTuple{<:Any,Int}}=nothing,
    kw...
) where K
    to_as_dims = _extent2dims(to; size, res, kw...)
    return Rasterizer(to_as_dims, data; kw..., fill)
end
function Rasterizer(to::DimTuple, data::T; fill, geomcolumn=nothing, kw...) where T
    if Tables.istable(T) # typeof so we dont check iterable table fallback in Tables.jl
        schema = Tables.schema(data)
        cols = Tables.columns(data)
        colnames = Tables.columnnames(Tables.columns(data))
        fillitr = _iterable_fill(cols, fill)
        # If fill is a symbol or tuple of Symbol we need to allocate based on the column type
        geomcolname = if isnothing(geomcolumn)
            geomcols = GI.geometrycolumns(data)
            isnothing(geomcols) ? nothing : first(geomcols)
        else
            geomcolumn
        end

        geometries = if geomcolname isa Symbol && geomcolname in Tables.columnnames(cols)
            # Its a geometry table
            Tables.getcolumn(cols, geomcolname)
        else
            # Its a point table
            pointcolnames = isnothing(geomcolumn) ? map(DD.dim2key, _auto_dim_columns(data, DEFAULT_POINT_ORDER)) : geomcolumn
            pointcols = map(k -> Tables.getcolumn(cols, k), pointcolnames)
            zip(pointcols...)
        end
        Rasterizer(to, geometries, fill, fillitr; kw...)
    else
        Rasterizer(to, GeoInterface.trait(data), data; fill, kw...)
    end
end
function Rasterizer(to::DimTuple, ::GI.AbstractFeatureCollectionTrait, fc; name, fill, kw...)
    # TODO: how to handle when there are fillvals with different types
    # fillval = _featurefillval(GI.getfeature(fc, 1), fill)
    fillitr = _iterable_fill(fc, fill)
    geometries = map(f -> GI.geometry(f), GI.getfeature(fc))
    Rasterizer(to, geometries; kw..., fill, fillitr)
end
function Rasterizer(to::DimTuple, ::GI.AbstractFeatureTrait, feature; fill, kw...)
    fillitr = _iterable_fill(feature, fill)
    # fillval = _featurefillval(feature, fill)
    Rasterizer(to, GI.geometry(feature), fill, fillitr; kw...)
end
function Rasterizer(to::DimTuple, ::Nothing, geoms; fill, kw...)
    fillitr = _iterable_fill(geoms, fill)
    Rasterizer(to, geoms, fill, fillitr; kw...)
end
function Rasterizer(to::DimTuple, ::GI.AbstractGeometryTrait, geom; fill, kw...)
    Rasterizer(to, geom, fill, _iterable_fill(geom, fill); kw...)
end

_fill_key_error(names, fill) = throw(ArgumentError("fill key $fill not found in table, use one of: $(Tuple(names))"))

# _featurefillval
# Get fill value from a feature, or use fill itself
_featurefillval(feature, fill::Nothing) = first(GI.properties(feature))
_featurefillval(feature, fill::Symbol) = GI.properties(feature)[fill]
_featurefillval(feature, fill::Val) = _featurefillval(feature, _unwrap(fill))
_featurefillval(feature, fill::NamedTuple) = map(f -> _featurefillval(feature, f), _unwrap(fill))
function _featurefillval(feature, fill::NTuple{<:Any,Symbol})
    map(fill) do key
        getproperty(GI.properties(feature), key)
    end |> NamedTuple{fill}
end
_featurefillval(feature, fill) = fill

_filter_name(name, fill::NamedTuple) = keys(fill)
_filter_name(name::NamedTuple, fill::NamedTuple) = keys(fill)
_filter_name(name::Nothing, fill::Nothing) = nothing
_filter_name(name::DimensionalData.NoName, fill::Union{Symbol,NTuple{<:Any,Symbol}}) = fill
_filter_name(name::Union{NamedTuple,Tuple,Array}, fill::NTuple{<:Any,Symbol}) = fill
function _filter_name(name::Union{NamedTuple,Tuple,Array}, fill::Union{Tuple,Array})
    length(name) == length(fill) || throw(ArgumentError("`name` keyword (possibly from `to` object) does not match length of fill. A fix is to use a `NamedTuple` for `fill`."))
    return name isa NamedTuple ? keys(name) : name
end
function _filter_name(name, fill)
    fill isa Union{Symbol,NTuple{<:Any,Symbol}} ? fill : name
end

# A Tuple of `Symbol` is multiple keys to make a RasterStack
_iterable_fill(data, keys::Tuple{Symbol,Vararg}) =
    NamedTuple{keys}(map(k -> _iterable_fill(data, k), keys))
# A Symbol is a Table or FeatureCollection key, it cant be used as fill itself
function _iterable_fill(data, key::Symbol)
    if GI.isfeature(data)
        return get(x -> error("feature has no key $key"), GI.properties(data), key)
    end
    cols = Tables.columns(data)
    # For column tables, get the column now
    names = Tables.columnnames(cols)
    key in names || _fill_key_error(names, key)
    return Tables.getcolumn(cols, key)
end
_iterable_fill(data, fill::Function) = fill
_iterable_fill(data, fill::NamedTuple) = begin
    map(f -> _iterable_fill(data, f), fill)
end
# Inspect our data and fill as much as possible to check they match
# and cycle any fill of known length one
function _iterable_fill(data, fill)
    if GI.isgeometry(data) || GI.isfeature(data) 
        return fill
    end
    if fill isa Number
        return Iterators.cycle(fill)
    end
    if Tables.istable(typeof(data))
        # we don't need the keys, just the column length
        data = first(Tables.columns(data))
    end
    if Base.IteratorSize(data) isa Union{Base.HasShape,Base.HasLength}
        fillvec = collect(fill)
        l = length(fillvec)
        n = length(data) 
        if l == 1
            # Cycle all length one iterables to fill every row
            return Iterators.cycle(fillvec[1])
        elseif !(l == n)
            throw(ArgumentError("Length of fill $l does not match length of iterator $n"))
        else
            return fillvec
        end
    else
        return fill
    end
end

_getfill(itrs::NamedTuple, i::Int) = map(itr -> _getfill(itr, i), itrs)
_getfill(itr::AbstractArray, i::Int) = itr[i]
_getfill(itr::Iterators.Cycle, i::Int) = first(itr)
_getfill(itr, i) = itr

"""
    rasterize([reduce], data; to, fill, kw...)

Rasterize a GeoInterface.jl compatable geometry or feature,
or a Tables.jl table with a `:geometry` column of GeoInterface.jl objects,
or `X`, `Y` points columns.

# Arguments

- `reduce`: a reducing function to reduce the fill value for all geometries that
    cover or touch a pixel down to a single value. The default is `last`.
    Any  that takes an iterable and returns a single value will work, including
    custom functions. However, there are optimisations for built-in methods
    including `sum`, `first`, `last`, `minimum`, `maximum`, `extrema` and `Statistics.mean`.
    These may be an order of magnitude or more faster than
    `count` is a special-cased as it does not need a fill value.
- `data`: a GeoInterface.jl `AbstractGeometry`, or a nested `Vector` of `AbstractGeometry`,
    or a Tables.jl compatible object containing a `:geometry` column or points and values columns.

# Keywords

These are detected automatically from `data` where possible.

$GEOM_KEYWORDS
- `fill`: the value or values to fill a polygon with. A `Symbol` or tuple of `Symbol` will
    be used to retrieve properties from features or column values from table rows. An array
    or other iterable will be used for each geometry, in order.
- `filename`: a filename to write to directly, useful for large files.
- `suffix`: a string or value to append to the filename.
    A tuple of `suffix` will be applied to stack layers. `keys(st)` are the default.
- `progress`: show a progress bar, `true` by default, `false` to hide..
- `geometrycolumn`: `Symbol` to manually select the column the geometries are in,
    or a tuple of `Symbol` for columns of point coordsates.

# Example

Rasterize a shapefile for China and plot, with a border.

```jldoctest
using Rasters, Plots, Dates, Shapefile, Downloads
using Rasters.LookupArrays

# Download a borders shapefile
shapefile_url = "https://github.com/nvkelso/natural-earth-vector/raw/master/10m_cultural/ne_10m_admin_0_countries.shp"
shapefile_name = "country_borders.shp"
isfile(shapefile_name) || Downloads.download(shapefile_url, shapefile_name)

# Load the shapes for china
china_border = Shapefile.Handle(shapefile_name).shapes[10]

# Rasterize the border polygon
china = rasterize(china_border; res=0.1, missingval=0, fill=1, boundary=:touches, progress=false)

# And plot
p = plot(china; color=:spring, legend=false)
plot!(p, china_border; fillalpha=0, linewidth=0.6)

savefig("build/china_rasterized.png"); nothing

# output

```

![rasterize](china_rasterized.png)

$EXPERIMENTAL
"""
function rasterize end
function rasterize(reduce::Function, data; kw...)
    rasterize(data; reduce, name=Symbol(string(reduce)), kw...)
end

_count_init_info(init) = @info "`rasterize` with `count` does not use the `init` keyword, $init ignored"
_count_fill_info(fill) = @info "`rasterize` with `count` does not use the `fill` keyword, $fill ignored"
_count_fill(x) = x + 1

# Catch some functions early

# count is faster with an incrementing function as `fill`
function rasterize(reduce::typeof(count), data; fill=nothing, init=nothing, kw...)
    isnothing(init) || _count_init_info(init)
    isnothing(fill) || _count_fill_info(fill)
    rasterize(data; kw..., name=:count, init=0, reduce=nothing, fill=_count_fill, missingval=0)
end
# `mean` is sum / count
# We can do better than this, but its easy for now
function rasterize(reduce::typeof(DD.Statistics.mean), data; fill, kw...)
    sums = rasterize(sum, data; kw..., fill)
    counts = rasterize(count, data; kw..., fill=nothing)
    rebuild(sums ./ counts; name=:mean)
end
function rasterize(data; to=nothing, fill, kw...)
    r = Rasterizer(to, data; fill, kw...)
    return Base.invokelatest() do
        return create_rasterize_dest(r) do dest
            _rasterize!(dest, r)
        end
    end
end

######################################
# Create a dest Array to rasterize! into

# create_rasterize_dest
# We create a Raster or RasterStack and apply f to it.
# This may be on disk, which is the reason for applying f rather than just
# returning the initiallised object - we may need to open it to be able to write.
create_rasterize_dest(f, r::Rasterizer) = create_rasterize_dest(f, r.eltype, r::Rasterizer)
# function _create_rasterize_dest(f, dims; fill, name=nothing, init=nothing, kw...)
    # _create_rasterize_dest(f, fill, init, name, dims; fill, kw...)
# end
function create_rasterize_dest(f, ::NamedTuple{K}, r::Rasterizer) where K
    layers = map(NamedTuple{K}(K), r.missingval, r.eltype) do name, missingval, eltype
        alloc_rasterize(r; eltype, name, missingval, suffix=name) do a
            # We should run `f` here, but it doesn't work yet.
            a
        end
    end
    # Combine layers into a RasterStack
    st = RasterStack(layers)
    # Apply f to the stack. (remove when we can do this in `_alloc_rasterize` while it is open)
    open(f, st; write=true)
    # Return the updated stack
    return st
end
function create_rasterize_dest(f, _, r::Rasterizer)
    result = alloc_rasterize(r) do a
        f(a)
    end
    return result
end

function alloc_rasterize(f, r::Rasterizer;
    eltype=r.eltype,
    name=r.name,
    missingval=r.missingval,
    metadata=r.metadata,
    suffix=r.suffix,
)
    A = create(r.filename, eltype, r.to; name, missingval, metadata, suffix)
    # TODO f should apply to the file when it is initially created
    # instead of reopening but we need a `create(f, filename, ...)` method
    open(A; write=true) do A
        A .= Ref(missingval)
        f(A)
    end
    return A
end

"""
    rasterize!(dest, data; fill)

Rasterize the geometries in `data` into the [`Raster`](@ref) or [`RasterStack`](@ref) `dest`,
using the values specified by `fill`.

# Arguments

- `dest`: a `Raster` or `RasterStack` to rasterize into.
- `data`: a GeoInterface.jl compatible object or an `AbstractVector` of such objects,
    or a Tables.jl compatible table containing `:geometry` column of GeoInterface compatible objects or
    columns with point names `X` and `Y`.
- `fill`: the value to fill a polygon with. A `Symbol` or tuple of `Symbol` will
    be used to retrieve properties from features or column values from table rows.

# Keywords

These are detected automatically from `A` and `data` where possible.

$GEOM_KEYWORDS
- `shape`: force `data` to be treated as `:polygon`, `:line` or `:point`, where possible
    Points can't be treated as lines or polygons, and lines may not work as polygons, but
    an attempt will be made.

And specifically for `shape=:polygon`:

- `boundary`: include pixels where the `:center` is inside the polygon, where
    the polygon `:touches` the pixel, or that are completely `:inside` the polygon.

# Example

```jldoctest
using Rasters, Plots, Dates, Shapefile, GeoInterface, Downloads
using Rasters.LookupArrays

# Download a borders shapefile
shapefile_url = "https://github.com/nvkelso/natural-earth-vector/raw/master/10m_cultural/ne_10m_admin_0_countries.shp"
shapefile_name = "country_borders.shp"
isfile(shapefile_name) || Downloads.download(shapefile_url, shapefile_name)

# Load the shapes for indonesia
indonesia_border = Shapefile.Handle(shapefile_name).shapes[1]

# Make an empty EPSG 4326 projected Raster of the area of Indonesia
dimz = X(Projected(90.0:0.1:145; sampling=Intervals(), crs=EPSG(4326))),
       Y(Projected(-15.0:0.1:10.9; sampling=Intervals(), crs=EPSG(4326)))

A = zeros(UInt32, dimz; missingval=UInt32(0))

# Rasterize each indonesian island with a different number. The islands are
# rings of a multi-polygon, so we use `GI.getring` to get them all separately.
islands = collect(GeoInterface.getring(indonesia_border))
rasterize!(A, islands; fill=1:length(islands), progress=false)

# And plot
p = plot(Rasters.trim(A); color=:spring)
plot!(p, indonesia_border; fillalpha=0, linewidth=0.7)

savefig("build/indonesia_rasterized.png"); nothing

# output

```

![rasterize](indonesia_rasterized.png)

$EXPERIMENTAL
"""
function rasterize! end
rasterize!(reduce::Function, x::RasterStackOrArray, data; kw...) =
    rasterize!(x::RasterStackOrArray, data; reduce, kw...)
function rasterize!(reduce::typeof(count), x::RasterStackOrArray, data; fill=nothing, init=nothing, kw...)
    isnothing(fill) || @info _count_fill_info(fill)
    isnothing(init) || @info _count_init_info(init)
    rasterize!(x::RasterStackOrArray, data; kw..., reduce=nothing, op=nothing, fill=_count_fill, init=0)
end
function rasterize!(x::RasterStackOrArray, data; kw...)
    r = Rasterizer(x, data; eltype=eltype(x), kw...)
    return _rasterize!(x, r)
end


function _rasterize!(x::RasterStackOrArray, r::Rasterizer)
    if r.shape == points
        _rasterize_points!(x, r)
    else
        _rasterize!(x, GI.trait(r.geom), r.geom, r.fillitr, r)
    end
    return x
end

# Single geometry to rasterize
function _rasterize!(A, ::GI.AbstractGeometryTrait, geom, fill, r::Rasterizer; allocs=r.allocs)
    (; op, init, missingval, lock, shape, boundary, verbose, progress) = r
    if r.shape === :point
        hasburned = _rasterize_points!(A, GI.trait(geom), geom, fill, r)
    else
        ext = _extent(geom)
        V = view(A, Touches(ext))
        length(V) > 0 || return false

        bools = _init_bools(V; metadata=metadata(A))
        boolmask!(bools, geom; allocs, lock, shape, boundary, verbose, progress)
        hasburned = any(bools)
        if hasburned
            # Avoid race conditions with a SectorLock
            isnothing(lock) || Base.lock(lock, V)
            _fill!(V, bools, fill, op, init, missingval)
            isnothing(lock) || Base.unlock(lock)
        end
    end
    return hasburned
end
# Fill points
function _rasterize!(A, trait::GI.AbstractPointTrait, point, fill, r::Rasterizer; allocs=nothing)
    # Avoid race conditions whern Point is in a mixed set of Geometries
    # isnothing(r.lock) || Base.lock(r.lock, A)
    hasburned = _fill_point!(A, trait, point; fill, r.lock)
    # isnothing(r.lock) || Base.unlock(r.lock)
    # for all points we avoid parallel rasterization completely - this method should not be hit often
    return hasburned
end
function _rasterize!(A, trait::Nothing, geoms, fill, r::Rasterizer; allocs=r.allocs)
    t1 = GI.trait(first(skipmissing(geoms)))
    if r.shape === :point
        return _rasterize_points!(A, trait, geoms, fill, r)
    else
        # Everything else is rasterized as line or polygon geometries
        return _rasterize_iterable!(_prepare_for_burning(A), geoms, r)
    end
end

# We rasterize all iterables from here
function _rasterize_iterable!(A, geoms, r::Rasterizer)
    (; reduce, op, fillitr) = r
    _rasterize_iterable!(A, geoms, reduce, op, fillitr, r)
end

function _rasterize_iterable!(
    A, geoms, reduce::Nothing, op::Nothing, fillitr::Iterators.Cycle, r::Rasterizer
)
    (; allocs, lock, shape, boundary, verbose, progress) = r
    # We dont need to iterate the fill, so just mask
    mask = boolmask(geoms; to=A, collapse=true, metadata=metadata(A), allocs, lock, shape, boundary, verbose, progress)
    # And broadcast the fill
    broadcast_dims!(A, A, mask) do v, m
        m ? fillitr.xs : v
    end
    return true
end
# Stacks
function _rasterize_iterable!(A, geoms, reduce::Nothing, op::Nothing,
    fillitr::NamedTuple{<:Any,Tuple{<:Iterators.Cycle,Vararg}}, r::Rasterizer
)
    (; allocs, shape, boundary, verbose, progress) = r
    # We dont need to iterate the fill, so just mask
    mask = boolmask(geoms; to=A, collapse=true, metadata=metadata(A), allocs, lock, shape, boundary, verbose, progress)
    foreach(A, fillitr) do A, f
        # And broadcast the fill
        broadcast_dims!(A, A, mask) do v, m
            m ? f.xs : v
        end
    end
    return true
end
# Simple iterator
function _rasterize_iterable!(A, geoms, reduce::Nothing, op::Nothing, fillitr, r::Rasterizer)
    range = _geomindices(geoms)
    burnchecks = _alloc_burnchecks(range)
    p = r.progress ? _progress(length(geoms); desc="Rasterizing...") : nothing
    # Threads.@threads
    for i in _geomindices(geoms)
        geom = _getgeom(geoms, i)
        ismissing(geom) && continue
        allocs = _get_alloc(r.allocs)
        fill = _getfill(fillitr, i)
        res = _rasterize!(A, GI.trait(geom), geom, fill, r; allocs)
        burnchecks[i] = res
        r.progress && ProgressMeter.next!(p)
    end
    _set_burnchecks(burnchecks, metadata(A), r.verbose)
    return any(burnchecks)
end

# Iterator with `reduce` or `op`
function _rasterize_iterable!(A, geoms, reduce, op::Nothing, fillitr, r::Rasterizer)
    # If there still isn't any `op`, reduce over all the values later rather than iteratively
    return _reduce_fill!(reduce, A, geoms, fillitr, r)
end
function _rasterize_iterable!(A, geoms, reduce, op, fillitr, r::Rasterizer)
    # Use op in a fast iterative reduction
    p = r.progress ? _progress(length(geoms); desc="Rasterizing...") : nothing
    range = _geomindices(geoms)
    burnchecks = _alloc_burnchecks(range)
    # Threads.@threads 
    for i in range
        geom = geoms[i]
        ismissing(geom) && continue
        allocs = _get_alloc(r.allocs)
        fill = _getfill(fillitr, i)
        burnchecks[i] = _rasterize!(A, GI.trait(geom), geom, fill, r; allocs)
        r.progress && ProgressMeter.next!(p)
    end
    _set_burnchecks(burnchecks, metadata(A), r.verbose)
    return any(burnchecks)
end

function _contains(extent::Extents.Extent, point)
    extent.X[1] <= point[1] < extent.X[2] && extent.Y[1] <= point[2] < extent.Y[2]
end


################################
# Fast point rasterization
#
# geoms is a iterator of points
_rasterize_points!(A, r::Rasterizer) =
    _rasterize_points!(A, GI.trait(r.geom), r.geom, r.fillitr, r)
function _rasterize_points!(A, ::GI.AbstractGeometryTrait, geom, fill, r::Rasterizer)
    points = GI.getpoint(geom)
    fill1 =_iterable_fill(points, fill)
    _rasterize_points!(A, nothing, points, fill1, r)
end
function _rasterize_points!(A, ::Nothing, geoms, fillitr, r::Rasterizer)
    (; reduce, op, missingval, init) = r
    t1 = GI.trait(first(skipmissing(geoms)))
    if !(t1 isa GI.PointTrait)
        # Recurse down until we hit points
        if r.fill isa Function
            for geom in _getgeom(geoms)
                _rasterize_points!(A, GI.trait(geom), geom, fillitr, r)
            end
        else
            if fillitr isa NamedTuple
                ntfill = _maybe_namedtuple_itr(fillitr)
                for (geom, f) in zip(_getgeom(geoms), ntfill)
                    _rasterize_points!(A, GI.trait(geom), geom, f, r)
                end
            else
                for (geom, f) in zip(_getgeom(geoms), fillitr)
                    _rasterize_points!(A, GI.trait(geom), geom, f, r)
                end
            end
        end
        return A
    end

    op = reduce == last ? _take_last : op
    if r.fill isa Function
        rasterize_points_unsorted!(A, geoms, nothing, fillitr, missingval, init)
    elseif op isa Function
        rasterize_points_unsorted!(A, geoms, op, fillitr, missingval, init)
    else
        rasterize_points_sorted!(A, geoms, reduce, fillitr, missingval, init)
    end
end
# Some algorithms don't need sort, like sum
function rasterize_points_unsorted!(A, points, op, fillitr, missingval, init)
    hasburned = false
    # Get extent information to properly shift the points
    # to the region of the array during rounding
    ext = Extents.extent(A)
    xrange = ext.X[2] - ext.X[1]
    yrange = ext.Y[2] - ext.Y[1]
    xsize = size(A, X)
    ysize = size(A, Y)
    n = 0
    for (i, point) in enumerate(points)
        ismissing(point) && continue
        point = (GI.x(point), GI.y(point))
        _contains(ext, point) || continue
        # Convert points to Int matching the array indices
        x = X(trunc(Int, (point[1] - ext.X[1]) / xrange * xsize) + 1)
        y = Y(trunc(Int, (point[2] - ext.Y[1]) / yrange * ysize) + 1)
        I = dims2indices(A, (x, y))
        _checkbounds(A, I...) || continue
        # Get the currenti fill value
        # And the Array value
        a = A[I...]
        f = _getfill(fillitr, i)
        if A isa RasterStack
            # Map for stacks
            f1 = map(a, f, init, missingval) do an, fn, it, mv
                _choose_fill(an, fn, op, it, mv)
            end
            A[I...] = f1
        else
            # Write directly for arrays
            f1 = _choose_fill(a, f, op, init, missingval)
            A[I...] = f1
        end
        # Mark that we have written at least one index
        hasburned = true
    end
    return hasburned
end

function rasterize_points_sorted!(A, points, reduce, fillitr, missingval, init)
    hasburned = false
    # Get extent information to properly shift the points
    # to the region of the array during rounding
    ext = Extents.extent(A)
    xrange = ext.X[2] - ext.X[1]
    yrange = ext.Y[2] - ext.Y[1]
    xsize = size(A, X)
    ysize = size(A, Y)
    # Convert all points to Int
    pf = map(points, _maybe_namedtuple_itr(fillitr)) do p, f
        ismissing(p) && return (0, 0) # Checkbounds will remove this later
        x = round(Int, (GI.x(p) - ext.X[1]) / xrange * xsize) + 1
        y = round(Int, (GI.y(p) - ext.Y[1]) / yrange * ysize) + 1
        ((x, y), f)
    end
    sort!(pf; by=first, alg=Base.Sort.DEFAULT_STABLE)
    prevpoint = first(first(pf))
    n = 0
    startind = 1
    for i in eachindex(pf)
        point, fill = pf[i]
        n += 1
        if prevpoint === point
            # We will reduce all of these points together later on
            continue
        else
            I = dims2indices(A, (X(prevpoint[1]), Y(prevpoint[2])))
            _checkbounds(A, I...) || continue
            startind = _fill_point!(reduce, A, I, pf, startind, i, missingval)
        end
        # Update the previous point to the current
        prevpoint = point
        # Mark that we have written at least one index
        hasburned = true
    end
    # Fill the last points
    I = dims2indices(A, (X(prevpoint[1]), Y(prevpoint[2])))
    _checkbounds(A, I...) && _fill_point!(reduce, A, I, pf, startind, lastindex(pf) + 1, missingval)
    return hasburned
end

# This will not be correct for multiple geometries
function _fill_point!(reduce, A, I, pf, startind, i, missingval)
    a = A[I...]
    v = _get_fill(reduce, pf[1][2], pf, startind:i - 1)
    x = _reduce_existing(reduce, a, v, missingval)
    A[I...] = x
    return i
end

_get_fill(reduce, ::NamedTuple{K}, pf, range) where K = NamedTuple{K}(map(k -> reduce(pf[n][2][k] for n in range), K))
_get_fill(reduce, ::Any, pf, range) = reduce(pf[n][2] for n in range)

_reduce_existing(reduce, as::NamedTuple, vs::NamedTuple, missingvals::NamedTuple) = 
    map((a, v, m) -> _reduce_existing(reduce, a, v, m), as, vs, missingvals)
function _reduce_existing(reduce, a, v, missingval)
    if ismissing(missingval) || (!ismissing(a) && a == missingval)
        v
    else
        reduce((a, v)) # This fill fail for mean, median etc
    end
end

# _reduce_fill!
#
# Mask `geoms` into each slice of a BitArray with the combined
# dimensions of `x` and `geoms`, then apply the supplied reducing
# function to each cell along the `:geometry` dimension. This uses a
# conditional iterator over the fill values and the corresponding
# rasterized Bool - true included the fill value in the iterator, false excludes it.
#
# We get 64 Bool values to a regular `Int` meaning this doesn't scale too
# badly for large tables of geometries. 64k geometries and a 1000 * 1000
# raster needs 1GB of memory just for the `BitArray`.
function _reduce_fill!(f, st::AbstractRasterStack, geoms, fill::NamedTuple, r::Rasterizer)
    (; allocs, lock, shape, boundary, verbose, progress) = r
    # Define mask dimensions, the same size as the spatial dims of x
    spatialdims = commondims(st, DEFAULT_POINT_ORDER)
    # Mask geoms as separate bool layers
    masks = boolmask(geoms; to=st, collapse=false, metadata=metadata(st), allocs, lock, shape, boundary, verbose, progress)
    # Use a generator over the array axis in case the iterator has no length
    geom_axis = axes(masks, Dim{:geometry}())
    fill = map(itr -> [v for (_, v) in zip(geom_axis, itr)], fill)
    T = NamedTuple{keys(st),Tuple{map(eltype, st)...}}
    _reduce_fill_inner!(f, st, T, geoms, fill, masks, progress)
end
function _reduce_fill!(f, A::AbstractRaster, geoms, fill, r::Rasterizer)
    (; allocs, lock, shape, boundary, verbose, progress) = r
    # Define mask dimensions, the same size as the spatial dims of x
    spatialdims = commondims(A, DEFAULT_POINT_ORDER)
    # Mask geoms as separate bool layers
    masks = boolmask(geoms; to=A, collapse=false, metadata=metadata(A), allocs, lock, shape, boundary, verbose, progress)
    # Use a generator over the array axis in case the iterator has no length
    geom_axis = parent(axes(masks, Dim{:geometry}()))
    fill = [val for (i, val) in zip(geom_axis, fill)]
    T = eltype(A)
    _reduce_fill_inner!(f, A, T, geoms, fill, masks, progress)
end

# Separated as a function barrier for type stability
function _reduce_fill_inner!(
    f, obj, ::Type{T}, geoms, fill::Union{AbstractVector,NamedTuple}, masks, progress::Bool
) where T
    p = progress ? _progress(size(obj, Y()); desc="Reducing...") : nothing
    Threads.@threads for y in axes(obj, Y())
        for x in axes(obj, X())
            D = (X(x), Y(y))
            # Do DimensionalData.jl indexing manually to avoid taking a view of the index and TwicePrecision problems
            I = dims2indices(masks, D)
            newval = _apply_reduction!(T, f, fill, view(parent(masks), I...))::Union{T,Nothing}
            if !isnothing(newval)
                @inbounds obj[D...] = newval::T
            end
        end
        progress && ProgressMeter.next!(p)
    end
    return obj
end

# _apply_reduction!
#
# Apply a reducing functin over an iterable
# This is applied for all reducing methods that don't have a matching `op` method
@inline function _apply_reduction!(::Type{T}, f, fill, pixel_geom_list) where T
    any(pixel_geom_list) || return nothing
    iterator = (fl for (fl, b) in zip(fill, pixel_geom_list) if b && !ismissing(fl))
    return convert(T, f(iterator))
end
@inline function _apply_reduction!(::Type{T}, f, fill::NamedTuple, pixel_geom_list) where T
    any(pixel_geom_list) || return nothing
    vals = map(fill) do fill
        iterator = (fl for (fl, b) in zip(fill, pixel_geom_list) if b && !ismissing(fl))
        f(iterator)
    end
    return convert(T, vals)
end


######################################
# Fill

_fill!(st::AbstractRasterStack, B, fill::NamedTuple, op, init::Union{Nothing,NamedTuple}, missingval) =
    _fill!(st, B, fill, op, init, map(_ -> missingval, fill))
function _fill!(st::AbstractRasterStack, B, fill::NamedTuple, op, init::Union{Nothing,NamedTuple}, missingval::NamedTuple)
    init = isnothing(init) ? map(_ -> nothing, fill) : init
    foreach(DimensionalData.layers(st), fill, init, missingval) do a, f, i, mv
        _fill!(a, B, f, op, i, mv)
    end
    return st
end
# If the array is initialised, we can use the existing values
function _fill!(A::AbstractRaster{T}, B, fill, op, init, missingval) where T
    broadcast_dims!(A, A, B) do a, b
        convert(T, b ? _choose_fill(a, fill, op, init, missingval) : a)::T
    end
    return A
end

@noinline function _choose_fill(::Type, a, b, fill::Function, op::Function, init, missingval)
    throw(ArgumentError("`fill` and `op` can't both be functions"))
end
@noinline function _choose_fill(::Type, a, b, fill::Function, op::Function, init, missingval::Missing)
    throw(ArgumentError("`fill` and `op` can't both be functions"))
end
# No op fill is a function, apply it unless missingval
_choose_fill(a, fill::Function, op::Nothing, init, missingval) =
    a == missingval ? fill(init) : fill(a)
# No op fill is a function, apply it unless missingval===missing
_choose_fill(a, fill::Function, op::Nothing, init, missingval::Missing) =
    ismissing(a) ? fill(init) : fill(a)
# No op fill is a function, no init
_choose_fill(a, fill::Function, op::Nothing, init::Nothing, missingval) = fill(a)
# No op fill is a function, no init fill a (repeated to avoid ambiguity)
_choose_fill(a, fill::Function, op::Nothing, init::Nothing, missingval::Missing) = fill(a)
# Op is a function, fill is not, missingval===missing
# apply retudcing op to a and fill, or to init and fill if a equals missing and init exists
function _choose_fill(a, fill, op, init, missingval::Missing)
    a1 = if ismissing(a)
        isnothing(init) ? a : init
    else
        a
    end
    _apply_op(op, a1, fill)
end
# Op is a function, fill is not, missingval===missing
# apply retudcing op to a and fill, or to init and fill if a equals missingval and init exists
function _choose_fill(a, fill, op, init, missingval)
    a1 = if a === missingval
        isnothing(init) ? a : init
    else
        a
    end
    _apply_op(op, a1, fill)
end

# apply reducing op to current value and fill value
_apply_op(op::Nothing, a1, fill) = fill
_apply_op(op::Function, a1, fill) = op(a1, fill)


_maybe_namedtuple_itr(nt::NamedTuple{K}) where K = (NamedTuple{K}(xs) for xs in zip(nt...))
_maybe_namedtuple_itr(itr) = itr
