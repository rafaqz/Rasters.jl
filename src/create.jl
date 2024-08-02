

"""
    create([filename], template::Raster; kw...)
    create([filename], T, template; kw...)

Create a new Raster. If `filename` is a `String` it will be created on disk,
and opened lazily. If it is `nothing` of not passed, a regular in-memory `Raster`
will be created. When written to disk, the values will be `missingval`,
if in-memory values will be `undef`.

The return value is a `Raster`. The `eltype` will usually be `T`, except
where `scale` and/or `offset` keywords are used, in which case `T` will
depend on the tyepe promotion of `scale` and `offset` with `T`.
`maskingval` will also affect the `eltype`.

## Arguments

- `filename`: a String file path, which will create a file on disk and return it as
    a lazy `Raster`, or `nothing` to create an in-memory `Raster`.
- `template`: a `Raster`, `Tuple` of `Dimension` or `Extents.Extent` to use as a template.
    If an `Extent` is used, a `size` or `res` keyword must be passed.
    If a `T` argument is not used, it is taken from the `template` eltype.
- `T`: the element type to use in the created array.

## Keywords

$NAME_KEYWORD
$REFDIMS_KEYWORD
$METADATA_KEYWORD
$WRITE_MISSINGVAL_KEYWORD
- `fillval`: A value to fill the array with. By default this will be
    `missingval`. If there is no `missingval` set or `fillval` is set to nothing
    disk values will remain undefined.
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
$CHUNKS_KEYWORD
- `reverse_y`: usually we want to write `Y` dimensions in reverse.
    When building dimensions from an `Extents.Extent` we do this by
    default, unless `reverse_y=false`. With template `Raster` or dimensions,
    the existing order is used.

## Example

Here we create a UInt8 GeoTIFF and open it as a Raster, from -80 to 80
lattitude, and 0 to 120 longitude, with a resolution of 0.25 degrees.

We scale values from 0-1 over `UInt8` 0-200, and using `255`.
Values that don't convert exactly will error (we could use `coerce=trunc` to fix that).

We use `UInt8(255)` as the `missingval` on disk, but mask it with `missing` in
the loaded `Raster`.

We use standard lat/lon (EPSG:4326) as the crs, and force writing if the file exists.

```julia
using Rasters, NCDatasets, ArchGDAL, Extents, Dates
using Rasters.Lookups
rast = Rasters.create("created.tif", UInt8, Extents.Extent(X=(0, 120), Y=(-80, 80), Band=(0, 12));
    res=(X=1.0, Y=1.0, Band=1),
    # size=(X=100, Y=100, Band=12),
    maskingval=nothing,
    name=:myraster,
    crs=EPSG(4326),
    force=true,
    sampling=(X=Intervals(Start()), Y=Intervals(Start()), Band=Intervals(Start())),
)
using ProfileView
@profview open(rast; write=true) do A
    A .= Rasters.Missings.nonmissingtype(eltype(A))(1)
    nothing
end
Raster("created.tif"; maskingval=nothing)
rm("created.tif")

extent = Extents.Extent(X=(0, 120), Y=(-80, 80))#, Band=(1, 3))
types = (a=UInt8, b=Int32, c=Float64=>Y)
rast = Rasters.create("created.nc", types, extent;
    # res=(X=1.0, Y=1.0, Band=1),
    maskingval=nothing,
    size=(X=100, Y=100),
    crs=EPSG(4326),
    force=true,
    sampling=(X=Intervals(Start()), Y=Intervals(Start()), Band=Points()),
)
RasterStack("created.nc")

╭───────────────────────────────────────────╮
│ 480×640 Raster{Union{Missing, Float64},2} │
├───────────────────────────────────────────┴───────────────────────────────────────── dims ┐
  ↓ X Projected{Float64} LinRange{Float64}(0.0, 119.75, 480) ForwardOrdered Regular Points,
  → Y Projected{Float64} LinRange{Float64}(79.75, -80.0, 640) ReverseOrdered Regular Points
├───────────────────────────────────────────────────────────────────────────────── metadata ┤
  Metadata{GDALsource} of Dict{String, Any} with 2 entries:
  "filepath" => "created.tif"
  "scale"    => 0.005
├─────────────────────────────────────────────────────────────────────────────────── raster ┤
  extent: Extent(X = (0.0, 119.75), Y = (-80.0, 79.75))
  missingval: missing
  crs: GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199
433,AUTHORITY["EPSG","9122"]],AXIS["Latitude",NORTH],AXIS["Longitude",EAST],AUTHORITY["EPSG","4326"]]
  filename: nothing
└───────────────────────────────────────────────────────────────────────────────────────────┘
```
"""
create(A::AbstractRaster; kw...) = create(nothing, A; kw...)
create(T::Union{Type,NamedTuple}, dims::Tuple; kw...) = create(nothing, T, dims; kw...)
create(T::Union{Type,NamedTuple}, extent::Extents.Extent; kw...) = create(nothing, T, dims; kw...)
create(filename::Union{AbstractString,Nothing}, A::AbstractRaster{T}; kw...) where T =
    create(filename, T, A; kw...)
function create(filename::Union{AbstractString,Nothing}, T::Union{Type,NamedTuple}, A::AbstractRaster;
    name=name(A),
    metadata=metadata(A),
    missingval=missingval(A),
    kw...
)
    return create(filename, T, dims(A); parent=parent(A), name, metadata, missingval, kw...)
end
function create(filename::AbstractString, T::Union{Type,NamedTuple}, dims::Tuple;
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
function create(filename::AbstractString, T::Union{Type,NamedTuple}, extent::Extents.Extent;
    res=nokw,
    size=nokw,
    crs=nothing,
    sampling=Points(),
    reverse_y=true,
    kw...
)
    ds = _extent2dims(extent; size, res, crs, sampling)
    ds = if reverse_y && hasdim(ds, Y())
        DD.setdims(ds, reverse(dims(ds, Y())))
    else
        ds
    end
    return create(filename, T, ds; kw...)
end
function create(filename::Nothing, ::Type{T}, dims::Tuple;
    parent=nokw,
    suffix=nokw,
    force=false,
    missingval,
    kw...
) where T
    eltype = isnothing(missingval) ? T : promote_type(T, typeof(missingval))
    data = if isnokw(parent) || isnothing(parent)
        Array{eltype}(undef, dims)
    else
        similar(parent, eltype, size(dims))
    end
    return Raster(data, dims; missingval, kw...)
end
function create(filename::Nothing, types::NamedTuple, dims::Tuple;
    suffix=nokw,
    force=false,
    missingval,
    kw...
)
    layers = map(types) do T
        # eltype = isnothing(missingval) ? T : promote_type(T, typeof(missingval))
        data = if isnokw(parent) || isnothing(parent)
            Array{eltype}(undef, dims)
        else
            similar(parent, eltype, size(dims))
        end
    end
    return RasterStack(layers, dims; missingval, kw...)
end
function create(filename::AbstractString, source::Source, ::Type{T}, dims::DimTuple;
    name=nokw,
    missingval=nokw,
    maskingval=missing,
    fillval=nokw,
    metadata=nokw,
    chunks=nokw,
    scale=nokw,
    offset=nokw,
    dropband=!hasdim(dims, Band),
    lazy=true,
    verbose=true,
    force=false,
    coerce=nokw,
) where T
    eltype = Missings.nonmissingtype(T)
    if isnokw(fillval) || isnothing(fillval)
        write = false # Leave fill undefined
        A = FillArrays.Zeros{eltype}(map(length, dims))
    else
        fillval isa T || throw(ArgumentError("fillval must be of type $T, got $fillval"))
        write = true # Write fill to disk
        A = FillArrays.Fill{eltype}(fillval, map(length, dims))
    end
    # Create layers of zero arrays
    rast = Raster(A, dims; name, missingval)
    Rasters.write(filename, source, rast;
        eltype, chunks, metadata, scale, offset, missingval, verbose, force, coerce, write
    )
    return Raster(filename; source, lazy, metadata, missingval, maskingval, dropband, coerce)
end
function create(filename::AbstractString, source::Source, layertypes::NamedTuple, dims::DimTuple;
    name=keys(layertypes),
    missingval=nokw,
    maskingval=missing,
    metadata=nokw,
    layerdims=nokw,
    layermetadata=nokw,
    chunks=nokw,
    scale=nokw,
    offset=nokw,
    dropband=!hasdim(dims, Band),
    lazy=true,
    verbose=true,
    force=false,
    coerce=nokw,
)
    layers = map(layertypes) do x
        if x isa Type
            eltype = Missings.nonmissingtype(x)
            size = map(length, dims)
        elseif x isa Pair{<:Type}
            eltype = Missings.nonmissingtype(x[1])
            ds = x[2]
            size = map(length, DD.dims(dims, DD._astuple(ds)))
        else
            throw(ArgumentError("Must be a Type or a Pair of Type and Dimension/Symbol"))
        end
        FillArrays.Zeros{eltype}(size)
    end
    layerdims = if isnokw(layerdims) 
        map(layertypes) do x
            if x isa Type
                DD.basedims(dims)
            else
                ds = DD._astuple(DD.basedims(x[2]))
            end
        end
    else
        layerdims
    end
    # if isnokw(fillval) || isnothing(fillval)
    #     write = false # Leave fill undefined
    #     A = FillArrays.Zeros{eltype}(map(length, dims))
    # else
    #     fillval isa T || throw(ArgumentError("fillval must be of type $T, got $fillval"))
    #     write = true # Write fill to disk
    #     A = FillArrays.Fill{eltype}(fillval, map(length, dims))
    # end
    # Create layers of zero arrays
    stack = RasterStack(layers, dims; layerdims, layermetadata, missingval)
    fn = Rasters.write(filename, stack;
        chunks, metadata, scale, offset, missingval, maskingval, verbose, force, coerce, write=false
    )
    return RasterStack(fn; source, lazy, metadata, layerdims, maskingval, dropband, coerce)
end
