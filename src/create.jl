

"""
    create([filename], template::Raster; kw...)
    create([filename], T, template::Raster; kw...)
    create([filename], T, template::Tuple; kw...)

Create a new Raster. If `filename` is a `String` it will be created on disk,
and opened lazily. If it is `nothing` a regular in-memory `Raster`
will be created. If written to disk, the values will be `missingval` when it
is defined, if in-memory values will be `undef`.

Generally all indices should be written to after `create`.

The return value is a `Raster`. The `eltype` will usually be `T`, except
where `scale` and/or `offset` keywords are used, in which case `T` will
depend on the tyepe promotion of `scale` and `offset` and `T`.
`maskingval` will also affect the `eltype`.

## Arguments

- `filename`: a String file path, which will create a file on disk and return it as
    a lazy `Raster`, or `nothing` to create an in-memory `Raster`.
- `T`: the element type to use in the created array.
- `template`: a `Raster`, `Tuple` of `Dimension` or `Extents.Extent` to use as a template. 
    If an `Extent` is used, a `size` or `res` keyword must be passed.

## Keywords

$NAME_KEYWORD
$REFDIMS_KEYWORD
$METADATA_KEYWORD
$MISSINGVAL_KEYWORD
$MASKINGVAL_KEYWORD
$SOURCE_KEYWORD
- `lazy`: A `Bool` specifying if to load data lazily from disk. For `create`
    `lazy=true` is the default, as creating a disk-based file is normally associated
    with it being larger than memory.
$CHUNKS_KEYWORD
$SCALE_KEYWORD
$OFFSET_KEYWORD
$COERCE_KEYWORD
$VERBOSE_KEYWORD
$RES_KEYWORD
$SIZE_KEYWORD
$CRS_KEYWORD
"""
create(A::AbstractRaster; kw...) where T = create(nothing, A; kw...)
create(T::Type, dims::Tuple; kw...) where T = create(nothing, T, dims; kw...)
create(T::Type, extent::Extents.Extent; kw...) where T = create(nothing, T, dims; kw...)
create(filename::Union{AbstractString,Nothing}, A::AbstractRaster{T}; kw...) where T =
    create(filename, T, A; kw...)
function create(filename::Union{AbstractString,Nothing}, T::Type, A::AbstractRaster;
    name=name(A),
    metadata=metadata(A),
    missingval=missingval(A),
    kw...
)
    return create(filename, T, dims(A); parent=parent(A), name, metadata, missingval, kw...)
end
function create(filename::AbstractString, T::Type, dims::Tuple;
    lazy=true,
    parent=nokw,
    suffix=nokw,
    source::Source=_sourcetrait(filename),
    missingval=nokw,
    kw...
)
    filename = _maybe_add_suffix(filename, suffix)
    # This calls `create` in the /sources file for this `source`
    return create(filename, source, T, dims; lazy, missingval, kw...)
end
function create(filename::Union{AbstractString,Extent}, T::Type, extent::Extents.Extent;
    res=nokw, 
    size=nokw, 
    crs=nothing,
    sampling=Points(),
    kw...
)
    ds = _extent2dims(extent; size, res, crs, sampling)
    return create(filename, T, ds; kw...)
end
function create(filename::Nothing, T::Type, dims::Tuple;
    parent=nokw,
    suffix=nokw,
    force=false,
    missingval,
    kw...
)
    T = isnothing(missingval) ? T : promote_type(T, typeof(missingval))
    data = isnokw(parent) || isnothing(parent) ? Array{T}(undef, dims) : similar(parent, T, size(dims))
    return Raster(data, dims; missingval, kw...)
end
function create(filename::AbstractString, source::Source, T::Type, dims::DimTuple;
    name=nokw,
    missingval=nokw,
    maskingval=missing,
    metadata=nokw,
    chunks=nokw,
    scale=nokw,
    offset=nokw,
    dropband=!hasdim(dims, Band),
    lazy=true,
    verbose=true,
    force=false,
    coerce=nokw,
)
    T1 = Missings.nonmissingtype(T)
    if isnothing(missingval)
        A = FillArrays.Zeros{T1}(map(length, dims))
    else
        missingval = ismissing(missingval) || isnokw(missingval) ? _type_missingval(T1) : convert(T1, missingval)
        A = FillArrays.Fill{T1}(missingval, map(length, dims))
    end
    # Create layers of zero arrays
    rast = Raster(A, dims; name, missingval)
    write(filename, source, rast; chunks, metadata, scale, offset, missingval, verbose, force, coerce)
    return Raster(filename; source, lazy, metadata, missingval, maskingval, dropband, coerce)
end
