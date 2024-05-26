

"""
    create(filename, A::Raster; kw...)
    create(filename, T, dims::Tuple; kw...)

Create a new Raster. If `filename` is a `String` it will be created on disk, 
and opened lazily. If it is `nothing` a regular in-memory `Raster` 
will be created. If written to disk, the values will be `missingval` when it
is defined, if in-memory values will be `undef`.

Generally all indices should be written to after `create`.

The return value is a `Raster`. The `eltype` will usually be `T`, except
where `scale` and/or `offset` keywords are used, in which case `T` will
depend on the tyepe promotion of `scale` and `offset` and `T`. 
`maskingval` will also affect the `eltype`.

# Keywords


"""
create(filename, A::AbstractRaster{T}; kw...) where T = create(filename, T, A; kw...)
function create(filename, x, A::AbstractRaster;
    name=name(A), 
    metadata=metadata(A), 
    missingval=missingval(A), 
    kw...
)
    create(filename, x, dims(A); parent=parent(A), name, metadata, missingval, kw...)
end
function create(filename::AbstractString, x, dims::Tuple;
    lazy=true,
    parent=nokw,
    suffix=nokw,
    source::Source=_sourcetrait(filename),
    missingval=nokw, 
    kw...
)
    filename = _maybe_add_suffix(filename, suffix)
    # This calls `create` in the /sources file for this `source`
    create(filename, source, x, dims; lazy, missingval, kw...)
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
    Raster(data, dims; missingval, kw...)
end
function create(filename::AbstractString, source::Source, T::Type, dims::DimTuple;
    name=nokw,
    missingval=nokw,
    maskingval=missingval,
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
    return Raster(filename; source, lazy, metadata, missingval, maskingval, dropband)
end
