export GRDarray, GRDstack, GRDdimMetadata, GRDarrayMetadata

const GRD_INDEX_ORDER = ForwardIndex()
const GRD_LON_ARRAY = ForwardArray()
const GRD_LAT_ARRAY = ReverseArray()
const GRD_BAND_ARRAY = ForwardArray()
const GRD_LON_RELATION = ForwardRelation()
const GRD_LAT_RELATION = ReverseRelation()
const GRD_BAND_RELATION= ForwardRelation()

# Metadata ########################################################################

"""
[`Metadata`](@ref) wrapper for `GRDarray` dimension metadata.
"""
struct GRDdimMetadata{K,V} <: DimMetadata{K,V}
    val::Dict{K,V}
end

"""
[`Metadata`](@ref) wrapper for `GRDarray` metadata.
"""
struct GRDarrayMetadata{K,V} <: ArrayMetadata{K,V}
    val::Dict{K,V}
end


# GRD attributes wrapper. Only used during file load, for dispatch.
struct GRDattrib{T,F,A}
    filename::F
    attrib::A
end
GRDattrib(filename::AbstractString) = begin
    filename = first(splitext(filename))
    lines = readlines(filename * ".grd")
    entries = filter!(x -> !isempty(x) && !(x[1] == '['), lines)
    attrib = Dict(Pair(string.(strip.(match(r"([^=]+)=(.*)", st).captures[1:2]))...) for st in entries)
    T = grd_datatype_translation[attrib["datatype"]]

    GRDattrib{T,typeof(filename),typeof(attrib)}(filename, attrib)
end

filename(grd::GRDattrib) = grd.filename
attrib(grd::GRDattrib) = grd.attrib

dims(grd::GRDattrib, crs=nothing, mappedcrs=nothing) = begin
    attrib = grd.attrib
    crs = crs isa Nothing ? ProjString(attrib["projection"]) : crs

    ncols, nrows, nbands = size(grd)

    xbounds = parse.(Float64, (attrib["xmin"], attrib["xmax"]))
    ybounds = parse.(Float64, (attrib["ymin"], attrib["ymax"]))

    xspan = (xbounds[2] - xbounds[1]) / ncols
    yspan = (ybounds[2] - ybounds[1]) / nrows

    # Not fully implemented yet
    latlon_metadata = GRDdimMetadata(Dict())

    latmode = Projected(
        order=Ordered(GRD_INDEX_ORDER, GRD_LAT_ARRAY, GRD_LAT_RELATION),
        span=Regular(yspan),
        sampling=Intervals(Start()),
        crs=crs,
        mappedcrs=mappedcrs,
    )
    lonmode = Projected(
        order=Ordered(GRD_INDEX_ORDER, GRD_LON_ARRAY, GRD_LON_RELATION),
        span=Regular(xspan),
        sampling=Intervals(Start()),
        crs=crs,
        mappedcrs=mappedcrs,
    )
    lat = Lat(LinRange(ybounds[1], ybounds[2] - yspan, nrows), latmode, latlon_metadata)
    lon = Lon(LinRange(xbounds[1], xbounds[2] - xspan, ncols), lonmode, latlon_metadata)
    band = Band(1:nbands; mode=Categorical(Ordered()))
    lon, lat, band
end

metadata(grd::GRDattrib, args...) = begin
    metadata = GRDarrayMetadata()
    for key in ("creator", "created", "history")
        val = get(grd.attrib, key, "")
        if val != ""
            metadata[key] = val
        end
    end
    metadata
end

missingval(grd::GRDattrib{T}) where T = begin
    mv = try
        parse(T, grd.attrib["nodatavalue"])
    catch
        @warn "No data $(missingval) is not convertible to data type $T. `missingval` set to NaN."
        missing
    end
end

name(grd::GRDattrib) = Symbol(get(grd.attrib, "layername", ""))


Base.eltype(::GRDattrib{T}) where T = T

Base.size(grd::GRDattrib) = begin
    ncols = parse(Int, grd.attrib["ncols"])
    nrows = parse(Int, grd.attrib["nrows"])
    nbands = parse(Int, grd.attrib["nbands"])
    # The order is backwards to jula array order
    ncols, nrows, nbands
end

Base.Array(grd::GRDattrib) = mmapgrd(Array, grd)


# Array ########################################################################

"""
    GRDarray(filename::String;
        mappedcrs=nothing,
        dims=(),
        refdims=(),
        name=nothing,
        missingval=nothing,
        metadata=nothing)

A [`DiskGeoArray`](@ref) that loads .grd files lazily from disk.

`GRDarray`s are always 3 dimensional, and have [`Lat`](@ref), [`Lon`](@ref) and
[`Band`](@ref) dimensions.

## Arguments

- `filename`: `String` pointing to a grd file. Extension is optional.

## Keyword Arguments

- `mappedcrs`: CRS format like `EPSG(4326)` used in `Selectors` like `Between` and `At`, and
  for plotting. Can be any CRS `GeoFormat` from GeoFormatTypes.jl, like `WellKnownText`.
- `name`: `String` name for the array, taken from the files `layername` attribute unless passed in.
- `dims`: `Tuple` of `Dimension`s for the array. Detected automatically, but can be passed in.
- `refdims`: `Tuple of` position `Dimension`s the array was sliced from.
- `missingval`: Value reprsenting missing values. Detected automatically when possible, but
  can be passed it.
- `metadata`: [`Metadata`](@ref) object for the array. Detected automatically as
  [`GRDarrayMetadata`](@ref), but can be passed in.

## Example

```julia
A = GRDarray("folder/file.grd"; mappedcrs=EPSG(4326))
# Select Australia using lat/lon coords, whatever the crs is underneath.
A[Lat(Between(-10, -43), Lon(Between(113, 153)))
```
"""
struct GRDarray{T,N,A,D<:Tuple,R<:Tuple,Na<:Symbol,Me,Mi,S
               } <: DiskGeoArray{T,N,D,LazyArray{T,N}}
    filename::A
    dims::D
    refdims::R
    name::Na
    metadata::Me
    missingval::Mi
    size::S
end
GRDarray(filename::String; kwargs...) =
    GRDarray(GRDattrib(filename), filename; kwargs...)
GRDarray(grd::GRDattrib, filename, key=nothing;
         crs=nothing,
         mappedcrs=nothing,
         dims=dims(grd, crs, mappedcrs),
         refdims=(),
         name=name(grd),
         missingval=missingval(grd),
         metadata=metadata(grd),
        ) = begin

    size_ = map(length, dims)

    T = eltype(grd)
    N = length(size_)
    name = Symbol(name)
    GRDarray{T,N,typeof.((filename,dims,refdims,name,metadata,missingval,size_))...
            }(filename, dims, refdims, name, metadata, missingval, size_)
end


# Base methods

"""
    Base.write(filename::AbstractString, ::Type{GRDarray}, s::AbstractGeoArray)

Write a [`GRDarray`](@ref) to a .grd file, with a .gri header file. The extension of
`filename` will be ignored.

Currently the `metadata` field is lost on `write` for `GRDarray`.

Returns `filename`.
"""
Base.write(filename::String, ::Type{<:GRDarray}, A::AbstractGeoArray) = begin
    if hasdim(A, Band)
        correctedA = permutedims(A, (Lon, Lat, Band)) |>
            a -> reorder(a, GRD_INDEX_ORDER) |>
            a -> reorder(a, (Lon(GRD_LON_RELATION), Lat(GRD_LAT_RELATION), Band(GRD_BAND_RELATION)))
        checkarrayorder(correctedA, (GRD_LON_ARRAY, GRD_LAT_ARRAY, GRD_BAND_ARRAY))
        nbands = length(val(dims(correctedA, Band)))
    else
        correctedA = permutedims(A, (Lon, Lat)) |>
            a -> reorder(a, GRD_INDEX_ORDER) |>
            a -> reorder(a, (Lon(GRD_LON_RELATION), Lat(GRD_LAT_RELATION)))
            checkarrayorder(correctedA, (GRD_LON_ARRAY, GRD_LAT_ARRAY))
        nbands = 1
    end
    # Remove extension
    filename = splitext(filename)[1]

    lon, lat = map(dims(A, (Lon(), Lat()))) do d
        convertmode(Projected, d)
    end
    ncols, nrows = size(A)
    xmin, xmax = bounds(lon)
    ymin, ymax = bounds(lat)
    proj = convert(String, convert(ProjString, crs(lon)))
    datatype = rev_grd_datatype_translation[eltype(A)]
    nodatavalue = missingval(A)
    minvalue = minimum(filter(x -> x !== missingval(A), data(A)))
    maxvalue = maximum(filter(x -> x !== missingval(A), data(A)))

    # Data: gri file
    open(filename * ".gri", "w") do IO
        write(IO, data(correctedA))
    end

    # Metadata: grd file
    open(filename * ".grd", "w") do IO
        write(IO,
            """
            [general]
            creator=GeoData.jl
            created= $(string(now()))
            [georeference]
            nrows= $nrows
            ncols= $ncols
            xmin= $xmin
            ymin= $ymin
            xmax= $xmax
            ymax= $ymax
            projection= $proj
            [data]
            datatype= $datatype
            nodatavalue= $nodatavalue
            byteorder= little
            nbands= $nbands
            minvalue= $minvalue
            maxvalue= $maxvalue
            [description]
            layername= $(name(A))
            """
        )
    end
    return filename
end


# AbstractGeoStack methods

"""
    GRDstack(filenames; keys, kwargs...)
    GRDstack(filenames...; keys, kwargs...)
    GRDstack(filenames::NamedTuple;
             window=(),
             metadata=nothing,
             childkwargs=(),
             refdims=())

Convenience method to create a DiskStack of [`GRDarray`](@ref) from `filenames`.

## Arguments

- `filenames`: A NamedTuple of stack keys and `String` filenames, or a `Tuple`,
  `Vector` or splatted arguments of `String` filenames.

## Keyword arguments

- `keys`: Used as stack keys when a `Tuple`, `Vector` or splat of filenames are passed in.
- `window`: A `Tuple` of `Dimension`/`Selector`/indices that will be applied to the
  contained arrays when they are accessed.
- `metadata`: Metadata as a [`StackMetadata`](@ref) object.
- `childkwargs`: A `NamedTuple` of keyword arguments to pass to the `childtype` constructor.
- `refdims`: `Tuple` of  position `Dimension` the array was sliced from.

## Example

Create a `GRDstack` from four files, that sets the child arrays
`mappedcrs` value when they are loaded.

```julia
files = (:temp="temp.tif", :pressure="pressure.tif", :relhum="relhum.tif")
stack = GRDstack(files; childkwargs=(mappedcrs=EPSG(4326),))
stack[:relhum][Lat(Contains(-37), Lon(Contains(144))
```
"""
GRDstack(args...; kwargs...) =
    DiskStack(args...; childtype=GRDarray, kwargs...)


withsource(f, ::Type{<:GRDarray}, filename::AbstractString, key...) =
    f(GRDattrib(filename))
withsourcedata(f, ::Type{<:GRDarray}, filename::AbstractString, key...) =
    mmapgrd(f, GRDattrib(filename))
withsourcedata(f, A::GRDarray, key...) =
    mmapgrd(f, A)

# Utils ########################################################################
mmapgrd(f, grd::Union{GRDarray,GRDattrib}) =
    mmapgrd(f, filename(grd), eltype(grd), size(grd))
mmapgrd(f, filename::AbstractString, T::Type, size::Tuple) = begin
    open(filename * ".gri", "r") do io
        mmap = Mmap.mmap(io, Array{T,length(size)}, size)
        output = f(mmap)
        close(io)
        output
    end
end

const grd_datatype_translation = Dict{String, DataType}(
    "LOG1S" => Bool,
    "INT1S" => Int8,
    "INT2S" => Int16,
    "INT4S" => Int32,
    "INT8S" => Int64,
    "INT1U" => UInt8,
    "INT2U" => UInt16,
    "FLT4S" => Float32,
    "FLT8S" => Float64
)

const rev_grd_datatype_translation =
    Dict{DataType, String}(v => k for (k,v) in grd_datatype_translation)
