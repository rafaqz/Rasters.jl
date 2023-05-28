_take_last(a, b) = b

_reduce_op(::typeof(sum)) = Base.add_sum
_reduce_op(::typeof(prod)) = Base.mul_prod
_reduce_op(::typeof(minimum)) = min
_reduce_op(::typeof(maximum)) = max
_reduce_op(x) = nothing

_reduce_init(reducer, st::AbstractRasterStack) = map(A -> _reduce_init(reducer, A), st)
_reduce_init(reducer, ::AbstractRaster{T}) where T = _reduce_init(reducer, T)
_reduce_init(reducer, nt::NamedTuple) = map(x -> _reduce_init(reducer, x), nt)
_reduce_init(f, x) = _reduce_init(f, typeof(x))

_reduce_init(::Nothing, x::Type{T}) where T = zero(T)
_reduce_init(f::Function, ::Type{T}) where T = zero(f((zero(nonmissingtype(T)), zero(nonmissingtype(T)))))
_reduce_init(::typeof(sum), ::Type{T}) where T = zero(nonmissingtype(T))
_reduce_init(::typeof(prod), ::Type{T}) where T = oneunit(nonmissingtype(T))
_reduce_init(::typeof(minimum), ::Type{T}) where T = typemax(nonmissingtype(T))
_reduce_init(::typeof(maximum), ::Type{T}) where T = typemin(nonmissingtype(T))

struct FillChooser{F,I,M}
    fill::F
    init::I
    missingval::M
end

struct RasterCreator{E,D,MD,MV,C,MC}
    eltype::E
    to::D
    filename::Union{String,Nothing}
    suffix::String
    name::Symbol
    metadata::MD
    missingval::MV
    crs::C
    mappedcrs::MC
end
function RasterCreator(to::DimTuple; 
    eltype,
    fill,
    missingval,
    filename=nothing,
    suffix="",
    res=nothing, # We shouldn't need this but coverage does
    crs=nothing,
    mappedcrs=nothing,
    name=nothing,
    metadata=Metadata(Dict()),
    kw...
)
    name = Symbol(_filter_name(name, fill))
    to = _as_intervals(to) # Only makes sense to rasterize to intervals
    RasterCreator(eltype, to, filename, suffix, name, metadata, missingval, crs, mappedcrs)
end
RasterCreator(to::AbstractRaster, data; kw...) = RasterCreator(dims(to); kw...)
RasterCreator(to::AbstractRasterStack, data; kw...) = RasterCreator(dims(to); name, kw...)
RasterCreator(to::Nothing, data; kw...) = RasterCreator(_extent(data); kw...)
RasterCreator(to, data; kw...) = RasterCreator(_extent(to); kw...)
function RasterCreator(to::Extents.Extent;
    res::Union{Nothing,Real,NTuple{<:Any,<:Real}}=nothing,
    size::Union{Nothing,Int,NTuple{<:Any,Int}}=nothing, kw...
)
    to_as_dims = _extent2dims(to; size, res, kw...)
    return RasterCreator(to_as_dims; kw...)
end


# Handles all data validations needed to
# run before rasterizing
struct Rasterizer{T,G,F,R,O,I,M}
    eltype::T
    geom::G
    fillitr::F
    reducer::R
    op::O
    init::I
    missingval::M
    lock::Union{SectorLocks,Nothing}
    shape::Symbol
    boundary::Symbol
    verbose::Bool
    progress::Bool
    threaded::Bool
end
function Rasterizer(geom, fill, fillitr;
    reducer=nothing,
    op=nothing,
    missingval=nothing,
    shape=nothing,
    eltype=nothing,
    init=nothing,
    boundary=:center,
    filename=nothing,
    verbose=true,
    progress=true,
    threaded=true,
    kw...
)
    # A single geometry does not need a reducing function 
    if !GI.isgeometry(geom)
        isnothing(reducer) && isnothing(op) && !(fill isa Function) && throw(ArgumentError("either reducer, op or fill must be a function"))
    end
 
    op = _reduce_op(reducer)

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
    init = isnothing(init) ? _reduce_init(reducer, filleltype) : init
    if shape == :points &&
        !GI.isgeometry(geom) &&
        !GI.trait(first(geom)) isa GI.PointTrait &&
        !(reducer in stable_reductions)
        @warn "currently `:points` rasterization of multiple non-`PointTrait` geometries may be innaccurate for `reducer` methods besides $stable_reductions. Make a Rasters.jl github issue if you need this to work"
    end
    eltype, missingval = get_eltype_missingval(eltype, missingval, fillitr, init, filename, op, reducer)
    lock = threaded ? SectorLocks() : nothing

    return Rasterizer(eltype, geom, fillitr, reducer, op, init, missingval, lock, shape, boundary, verbose, progress, threaded)
end
function Rasterizer(data::T; fill, geomcolumn=nothing, kw...) where T
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
        Rasterizer(geometries, fill, fillitr; kw...)
    else
        Rasterizer(GeoInterface.trait(data), data; fill, kw...)
    end
end
function Rasterizer(::GI.AbstractFeatureCollectionTrait, fc; name, fill, kw...)
    # TODO: how to handle when there are fillvals with different types
    # fillval = _featurefillval(GI.getfeature(fc, 1), fill)
    fillitr = _iterable_fill(fc, fill)
    geometries = map(f -> GI.geometry(f), GI.getfeature(fc))
    Rasterizer(geometries; kw..., fill, fillitr)
end
function Rasterizer(::GI.AbstractFeatureTrait, feature; fill, kw...)
    fillitr = _iterable_fill(feature, fill)
    # fillval = _featurefillval(feature, fill)
    Rasterizer(GI.geometry(feature), fill, fillitr; kw...)
end
function Rasterizer(::Nothing, geoms; fill, kw...)
    fillitr = _iterable_fill(geoms, fill)
    Rasterizer(geoms, fill, fillitr; kw...)
end
function Rasterizer(::GI.AbstractGeometryTrait, geom; fill, kw...)
    Rasterizer(geom, fill, _iterable_fill(geom, fill); kw...)
end

function get_eltype_missingval(eltype, missingval, fillitr, init::NamedTuple, filename, op, reducer)
    eltype = eltype isa NamedTuple ? eltype : map(_ -> eltype, init)
    missingval = missingval isa NamedTuple ? missingval : map(_ -> missingval, init)
    em = map(eltype, missingval, fillitr, init) do et, mv, fi, i
        get_eltype_missingval(et, mv, fi, i, filename, op, reducer)
    end
    eltype = map(first, em)
    missingval = map(last, em)
    return eltype, missingval
end
function get_eltype_missingval(known_eltype, missingval, fillitr, init, filename, op, reducer)
    eltype = if isnothing(known_eltype)
        if fillitr isa Function
            typeof(fillitr(init))
        elseif op isa Function
            typeof(op(init, init))
        elseif reducer isa Function
            typeof(reducer((init, init)))
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


const RASTERIZE_KEYWORDS = """
- `fill`: the value or values to fill a polygon with. A `Symbol` or tuple of `Symbol` will
    be used to retrieve properties from features or column values from table rows. An array
    or other iterable will be used for each geometry, in order. `fill` can also be a function of 
    the current value, e.g. `x -> x + 1`.
- `op`: A reducing function that accepts two values and returns one, like `min` to `minimum`.
    For common methods this will be assigned for you, or is not required. But you can use it
    instead of a `reducer` as it will usually be faster.
- `shape`: force `data` to be treated as `:polygon`, `:line` or `:point`, where possible
    Points can't be treated as lines or polygons, and lines may not work as polygons, but
    an attempt will be made.
- `geometrycolumn`: `Symbol` to manually select the column the geometries are in
    when `data` is a Tables.jl compatible table, or a tuple of `Symbol` for columns of
    point coordinates.
- `progress`: show a progress bar, `true` by default, `false` to hide..
- `verbose`: print information and warnings whne there are problems with the rasterisation.
    `true` by default.
- `threaded`: run operations in parallel. `true` by default.
"""

const RASTERIZE_ARGUMENTS = """
- `reducer`: a reducing function to reduce the fill value for all geometries that
    cover or touch a pixel down to a single value. The default is `last`.
    Any  that takes an iterable and returns a single value will work, including
    custom functions. However, there are optimisations for built-in methods
    including `sum`, `first`, `last`, `minimum`, `maximum`, `extrema` and `Statistics.mean`.
    These may be an order of magnitude or more faster than
    `count` is a special-cased as it does not need a fill value.
- `data`: a GeoInterface.jl `AbstractGeometry`, or a nested `Vector` of `AbstractGeometry`,
    or a Tables.jl compatible object containing a `:geometry` column or points and values columns.
"""

"""
    rasterize([reducer], data; kw...)

Rasterize a GeoInterface.jl compatable geometry or feature,
or a Tables.jl table with a `:geometry` column of GeoInterface.jl objects,
or `X`, `Y` points columns.

# Arguments

$RASTERIZE_ARGUMENTS

# Keywords

These are detected automatically from `data` where possible.

$GEOM_KEYWORDS
$RASTERIZE_KEYWORDS
- `filename`: a filename to write to directly, useful for large files.
- `suffix`: a string or value to append to the filename.
    A tuple of `suffix` will be applied to stack layers. `keys(st)` are the default.

# Example

Rasterize a shapefile for China and plot, with a border.

```jldoctest
using Rasters, RasterDataSources, ArchGDAL, Plots, Dates, Shapefile, Downloads
using Rasters.LookupArrays

# Download a borders shapefile
shapefile_url = "https://github.com/nvkelso/natural-earth-vector/raw/master/10m_cultural/ne_10m_admin_0_countries.shp"
shapefile_name = "country_borders.shp"
isfile(shapefile_name) || Downloads.download(shapefile_url, shapefile_name)

# Load the shapes for china
china_border = Shapefile.Handle(shapefile_name).shapes[10]

# Rasterize the border polygon
china = rasterize(last, china_border; res=0.1, missingval=0, fill=1, boundary=:touches, progress=false)

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
function rasterize(reducer::Function, data; kw...)
    rasterize(data; reducer, name=Symbol(string(reducer)), kw...)
end

_count_init_info(init) = @info "`rasterize` with `count` does not use the `init` keyword, $init ignored"
_count_fill_info(fill) = @info "`rasterize` with `count` does not use the `fill` keyword, $fill ignored"
_count_fill(x) = x + 1

# Catch some functions early

# count is faster with an incrementing function as `fill`
function rasterize(reducer::typeof(count), data; fill=nothing, init=nothing, kw...)
    isnothing(init) || _count_init_info(init)
    isnothing(fill) || _count_fill_info(fill)
    rasterize(data; kw..., name=:count, init=0, reducer=nothing, fill=_count_fill, missingval=0)
end
# `mean` is sum / count
# We can do better than this, but its easy for now
function rasterize(reducer::typeof(DD.Statistics.mean), data; fill, kw...)
    sums = rasterize(sum, data; kw..., fill)
    counts = rasterize(count, data; kw..., fill=nothing)
    rebuild(sums ./ counts; name=:mean)
end
function rasterize(data; to=nothing, fill, threaded=true, kw...)
    r = Rasterizer(data; fill, threaded, kw...)
    rc = RasterCreator(to, data; kw..., eltype=r.eltype, fill, missingval=r.missingval)
    allocs = r.shape == :points ? nothing : _burning_allocs(rc.to; threaded)
    return create_rasterize_dest(rc) do dest
        _rasterize!(dest, r; allocs)
    end
end

######################################
# Create a dest Array to rasterize! into

# create_rasterize_dest
# We create a Raster or RasterStack and apply f to it.
# This may be on disk, which is the reason for applying f rather than just
# returning the initiallised object - we may need to open it to be able to write.
create_rasterize_dest(f, r::RasterCreator) = create_rasterize_dest(f, r.eltype, r)
# function _create_rasterize_dest(f, dims; fill, name=nothing, init=nothing, kw...)
    # _create_rasterize_dest(f, fill, init, name, dims; fill, kw...)
# end
function create_rasterize_dest(f, ::NamedTuple{K}, r::RasterCreator) where K
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
function create_rasterize_dest(f, _, r::RasterCreator)
    result = alloc_rasterize(r) do a
        f(a)
    end
    return result
end

function alloc_rasterize(f, r::RasterCreator;
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
    rasterize!([reducer], dest, data; kw...)

Rasterize the geometries in `data` into the [`Raster`](@ref) or [`RasterStack`](@ref) `dest`,
using the values specified by `fill`.

# Arguments

- `dest`: a `Raster` or `RasterStack` to rasterize into.
$RASTERIZE_ARGUMENTS

# Keywords

These are detected automatically from `A` and `data` where possible.

$RASTERIZE_KEYWORDS
$GEOM_KEYWORDS

# Example

```jldoctest
using Rasters, RasterDataSources, ArchGDAL, Plots, Dates, Shapefile, GeoInterface, Downloads
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
rasterize!(last, A, islands; fill=1:length(islands), progress=false)

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
rasterize!(reducer::Function, x::RasterStackOrArray, data; kw...) =
    rasterize!(x::RasterStackOrArray, data; reducer, kw...)
function rasterize!(reducer::typeof(count), x::RasterStackOrArray, data; fill=nothing, init=nothing, kw...)
    isnothing(fill) || @info _count_fill_info(fill)
    isnothing(init) || @info _count_init_info(init)
    rasterize!(x::RasterStackOrArray, data; kw..., reducer=nothing, op=nothing, fill=_count_fill, init=0)
end
function rasterize!(x::RasterStackOrArray, data; threaded=true, kw...)
    r = Rasterizer(data; eltype=eltype(x), threaded, kw...)
    allocs = r.shape == :points ? nothing : _burning_allocs(dims(x); threaded)
    return _rasterize!(x, r; allocs)
end
function _rasterize!(A::RasterStackOrArray, r::Rasterizer; allocs=nothing)
    A1 = _prepare_for_burning(A)
    if r.shape == points
        _rasterize_points!(A1, r; allocs)
    else
        _rasterize!(A1, GI.trait(r.geom), r.geom, r.fillitr, r; allocs)
    end
    return A
end
# Single geometry to rasterize
function _rasterize!(A, ::GI.AbstractGeometryTrait, geom, fill, r::Rasterizer; allocs=nothing)
    (; op, init, missingval, lock, shape, boundary, verbose, progress) = r
    if r.shape === :point
        hasburned = _rasterize_points!(A, GI.trait(geom), geom, fill, r)
    else
        ext = _extent(geom)
        V = view(A, Touches(ext))
        length(V) > 0 || return false

        bools = _init_bools(commondims(V, DEFAULT_POINT_ORDER), Bool; metadata=metadata(A))
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
function _rasterize!(A, trait::Nothing, geoms, fill, r::Rasterizer; allocs=nothing)
    if r.shape === :point
        return _rasterize_points!(A, geoms, fill, r)
    else
        (; reducer, op, fillitr) = r
        # Everything else is rasterized as line or polygon geometries
        return _rasterize_iterable!(A, geoms, reducer, op, fillitr, r, allocs)
    end
end

# We rasterize all iterables from here
function _rasterize_iterable!(A, geoms, reducer, op, fillitr, r::Rasterizer, allocs)
    # reduce 
    if isnothing(op) && !(fillitr isa Function)
        return _reduce_bitarray!(reducer, A, geoms, fillitr, r, allocs)
    end
    range = _geomindices(geoms)
    burnchecks = _alloc_burnchecks(range)
    _run(range, r.threaded, r.progress, "Rasterizing...") do i
        geom = geoms[i]
        ismissing(geom) && return nothing
        a = _get_alloc(allocs)
        fill = _getfill(fillitr, i)
        burnchecks[i] = _rasterize!(A, GI.trait(geom), geom, fill, r; allocs=a)
        return nothing
    end
    _set_burnchecks(burnchecks, metadata(A), r.verbose)
    return any(burnchecks)
end

################################
# Fast point rasterization
#
# geoms is a iterator of points
_rasterize_points!(A, r::Rasterizer) = _rasterize_points!(A, r.geom, r.fillitr, r)
_rasterize_points!(A, geom, fillitr, r::Rasterizer) =
    _rasterize_points!(A, GI.trait(geom), geom, fillitr, r)
function _rasterize_points!(A, ::GI.AbstractGeometryTrait, geom, fill, r::Rasterizer)
    points = GI.getpoint(geom)
    fill1 =_iterable_fill(points, fill)
    _rasterize_points!(A, nothing, points, fill1, r)
end
function _rasterize_points!(A, ::Nothing, geoms, fillitr, r::Rasterizer)
    (; reducer, op, missingval, init) = r
    t1 = GI.trait(first(skipmissing(geoms)))
    hasburned = false
    if !(t1 isa GI.PointTrait)
        # Recurse down until we hit points
        if r.fillitr isa Function
            for geom in _getgeom(geoms)
                hasburned |= _rasterize_points!(A, GI.trait(geom), geom, fillitr, r)
            end
        else
            if fillitr isa NamedTuple
                ntfill = _maybe_namedtuple_itr(fillitr)
                for (geom, f) in zip(_getgeom(geoms), ntfill)
                    hasburned |= _rasterize_points!(A, GI.trait(geom), geom, f, r)
                end
            else
                for (geom, f) in zip(_getgeom(geoms), fillitr)
                    hasburned |= _rasterize_points!(A, GI.trait(geom), geom, f, r)
                end
            end
        end
        return hasburned
    end
    # Get extent information to properly shift the points
    # to the region of the array during rounding
    ext = Extents.extent(A)
    xrange = ext.X[2] - ext.X[1]
    yrange = ext.Y[2] - ext.Y[1]
    xsize = size(A, X)
    ysize = size(A, Y)
    s = (; ext, xrange, yrange, xsize, ysize)
    _rasterize_points_inner!(A, geoms, fillitr, s, reducer, op, missingval, init)
end

@noinline function _rasterize_points_inner!(A, geoms, fillitr::F, s, reducer::R, op::O, missingval, init)::Bool where {F,O,R}
    function xy(p) 
        # TODO handle reversed lookups
        x = round(Int, (GI.x(p) - s.ext.X[1]) / s.xrange * s.xsize) + 1
        y = round(Int, (GI.y(p) - s.ext.Y[1]) / s.yrange * s.ysize) + 1
        (x, y)
    end

    op = reducer == last ? _take_last : op
    if fillitr isa Function
        # We don't need to iterate fill
        points = Iterators.map(xy, geoms)
        return rasterize_points_fillfunc!(fillitr, A, points, missingval, init)
    end

    if op isa Function
        # `Iterators.map` allocates less
        points_fill = Iterators.map(geoms, _maybe_namedtuple_itr(fillitr)) do p, f
            (xy(p), f)
        end
        return rasterize_points_op!(op, A, points_fill, missingval, init)
    else
        # We need to use regular `map` to a vector for sorting later
        points_fill = map(geoms, _maybe_namedtuple_itr(fillitr)) do p, f
            (xy(p), f)
        end
        return rasterize_points_reduce!(reducer, A, points_fill, missingval, init)
    end
end

# Some algorithms don't need sort, like sum
@noinline function rasterize_points_fillfunc!(fillfunc::F, A, points, missingval, init)::Bool where F<:Function
    hasburned = false
    n = 0
    for point in points
        I = dims2indices(A, (X(point[1]), Y(point[2])))
        _checkbounds(A, I...) || continue
         _fill_func!(fillfunc, A, I)
        # Mark that we have written at least one index
        hasburned = true
    end
    return hasburned
end

function _fill_func!(fillfunc, A::Raster, I)
    @inbounds A[I...] = fillfunc(A[I...])
end
function _fill_func!(fillfunc, A::RasterStack, I)
    @inbounds a = A[I...]
    f1 = map(a) do x
        fillfunc(x)
    end
    @inbounds A[I...] = f1
end

# Some reductions don't need sort, like sum
@noinline function rasterize_points_op!(op, A, points_fill, missingval, init)::Bool
    hasburned = false
    n = 0
    
    for (point, fill) in points_fill
        I = dims2indices(A, (X(point[1]), Y(point[2])))
        _checkbounds(A, I...) || continue
         _fill_op!(op, A, fill, init, missingval, I)
        # Mark that we have written at least one index
        hasburned = true
    end
    return hasburned
end

function _fill_op!(op::O, A::Raster, fill, init, missingval, I) where {O<:Function}
    @inbounds a = A[I...]
    f1 = _choose_fill(op, a, FillChooser(fill, init, missingval))
    @inbounds A[I...] = f1
end
function _fill_op!(op::O, A::RasterStack, fill, init, missingval, I) where {O<:Function}
    @inbounds a = A[I...]
    choosers = map(FillChooser, fill, init, missingval)
    f1 = map(a, choosers) do an, fc
        _choose_fill(op, an, fc)
    end
    @inbounds A[I...] = f1
end

type_length(tup::Type{T}) where {T<:Union{Tuple,NamedTuple}} = length(tup.types)

@noinline function rasterize_points_reduce!(reducer, A, points_fill, missingval, init)
    hasburned = false
    # Convert all points to Int
    sort!(points_fill; by=first, alg=Base.Sort.DEFAULT_STABLE)
    prevpoint = first(first(points_fill))
    startind = 1
    for n in eachindex(points_fill)
        point, fill = points_fill[n]
        if prevpoint === point
            continue # We will reduce all of these points together later on
        else
            I = dims2indices(A, (X(prevpoint[1]), Y(prevpoint[2])))
            _checkbounds(A, I...) || continue
            startind = _fill_reduce!(reducer, A, I, points_fill, startind, n, missingval)
        end
        prevpoint = point # Update the previous point to the current
        hasburned = true # Mark that we have written at least one index
    end
    # Fill the last points
    I = dims2indices(A, (X(prevpoint[1]), Y(prevpoint[2])))
    n = lastindex(points_fill) + 1
    _checkbounds(A, I...) && _fill_reduce!(reducer, A, I, points_fill, startind, n, missingval)
    return hasburned
end

# This will not be correct for multiple geometries
function _fill_reduce!(reducer, A, I, pf, startind, i, missingval)
    @inbounds a = A[I...]
    v = _get_fill(reducer, pf[1][2], pf, startind:i - 1)
    x = _reduce_existing(reducer, a, v, missingval)
    @inbounds A[I...] = x
    return i
end

_get_fill(reducer, ::NamedTuple{K}, pf, range) where K = NamedTuple{K}(map(k -> reducer(pf[n][2][k] for n in range), K))
_get_fill(reducer, ::Any, pf, range) = reducer(pf[n][2] for n in range)

_reduce_existing(reducer, as::NamedTuple, vs::NamedTuple, missingvals::NamedTuple) =
    map((a, v, m) -> _reduce_existing(reducer, a, v, m), as, vs, missingvals)
function _reduce_existing(reducer, a, v, missingval)
    if ismissing(missingval) || (!ismissing(a) && a == missingval)
        v
    else
        reducer((a, v)) # This fill fail for mean, median etc
    end
end

# _reduce_bitarray!
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
function _reduce_bitarray!(f, st::AbstractRasterStack, geoms, fill::NamedTuple, r::Rasterizer, allocs)
    (; lock, shape, boundary, verbose, progress, threaded) = r
    # Define mask dimensions, the same size as the spatial dims of x
    spatialdims = commondims(st, DEFAULT_POINT_ORDER)
    # Mask geoms as separate bool layers
    masks = boolmask(geoms; to=st, collapse=false, metadata=metadata(st), allocs, lock, shape, boundary, verbose, progress)
    # Use a generator over the array axis in case the iterator has no length
    geom_axis = axes(masks, Dim{:geometry}())
    fill = map(itr -> [v for (_, v) in zip(geom_axis, itr)], fill)
    T = NamedTuple{keys(st),Tuple{map(eltype, st)...}}
    range = axes(first(st), Y())
    _run(range, threaded, progress, "Reducing...") do y
        _reduce_bitarray_loop(f, st, T, fill, masks, y)
    end
end
function _reduce_bitarray!(f, A::AbstractRaster, geoms, fill, r::Rasterizer, allocs)
    (; lock, shape, boundary, verbose, progress, threaded) = r
    # Define mask dimensions, the same size as the spatial dims of x
    spatialdims = commondims(A, DEFAULT_POINT_ORDER)
    # Mask geoms as separate bool layers
    masks = boolmask(geoms; to=A, collapse=false, metadata=metadata(A), lock, shape, boundary, verbose, progress, threaded)
    # Use a generator over the array axis in case the iterator has no length
    geom_axis = parent(axes(masks, Dim{:geometry}()))
    fill = [val for (i, val) in zip(geom_axis, fill)]
    T = eltype(A)
    range = axes(A, Y())
    _run(range, threaded, progress, "Reducing...") do y
        _reduce_bitarray_loop(f, A, T, fill, masks, y)
    end
    return A
end

function _reduce_bitarray_loop(f, A, ::Type{T}, fill, masks, y) where T
    for x in axes(A, X())
        D = (X(x), Y(y))
        # Do DimensionalData.jl indexing manually to avoid taking a view of the index and TwicePrecision problems
        I = dims2indices(masks, D)
        newval = _apply_reduction!(T, f, fill, view(parent(masks), I...))::Union{T,Nothing}
        if !isnothing(newval)
            @inbounds A[D...] = newval::T
        end
    end
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
        convert(T, b ? _choose_fill(op, a, FillChooser(fill, init, missingval)) : a)::T
    end
    return A
end

# No op fill is a function, apply it unless missingval
Base.@assume_effects :total _choose_fill(op::Nothing, a, fc::FillChooser{<:Function}) =
    a == missingval ? fc.fill(fc.init) : fc.fill(a)
# No op fill is a function, apply it unless missing
Base.@assume_effects :total _choose_fill(op::Nothing, a, fc::FillChooser{<:Function,<:Any,Missing}) =
    ismissing(a) ? fc.fill(fc.init) : fc.fill(a)
# No op fill is a function, no init
Base.@assume_effects :total _choose_fill(op::Nothing, a, fc::FillChooser{<:Function,Nothing}) = fc.fill(a)
# No op fill is a function, no init fill a (repeated to avoid ambiguity)
Base.@assume_effects :total _choose_fill(op::Nothing, a, fc::FillChooser{<:Function,Nothing,Missing}) = fc.fill(a)
# Op is a function, fill is not, missingval===missing
# apply retudcing op to a and fill, or to init and fill if a equals missing and init exists
Base.@assume_effects :total function _choose_fill(op::F, a, fc::FillChooser{<:Any,<:Any,Missing}) where F<:Function
    a1 = ismissing(a) ? fc.init : a
    _apply_op(op, a1, fc.fill)
end
Base.@assume_effects :total function _choose_fill(op::F, a, fc::FillChooser) where F<:Function
    a1 = a === fc.missingval ? fc.init : a
    _apply_op(op, a1, fc.fill)
end
Base.@assume_effects :total function _choose_fill(op::F, a, fc::FillChooser{<:Any,Nothing,Missing}) where F<:Function
    _apply_op(op, a, fc.fill)
end
# No op, fill is a value - this is just one geometry
Base.@assume_effects :total function _choose_fill(op, a, fc::FillChooser)
    fc.fill
end
# Op is a function, fill is not, missingval===missing
# apply retudcing op to a and fill, or to init and fill if a equals missingval and init exists
# @inline function _choose_fill(a, fill, op::F, init::Nothing, missingval) where F<:Function
#     _apply_op(op, a, fill)
# end

# apply reducing op to current value and fill value
Base.@assume_effects :total _apply_op(op::Nothing, a1, fill) = fill
Base.@assume_effects :total _apply_op(op::F, a1, fill) where F<:Function = op(a1, fill)


_maybe_namedtuple_itr(nt::NamedTuple{K}) where K = 
    (NamedTuple{K}(xs) for xs in zip(nt...)) 
_maybe_namedtuple_itr(itr) = itr
