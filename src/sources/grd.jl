export GrdArray, GrdStack, GrdMetadata, GrdDimMetadata

# Metadata ########################################################################

"""
[`ArrayMetadata`](@ref) wrapper for `GrdArray`.
"""
struct GrdMetadata{K,V} <: ArrayMetadata{K,V}
    val::Dict{K,V}
end
GrdMetadata() = GrdMetadata(Dict())

"""
[`DimMetadata`](@ref) wrapper for `GrdArray` dimensions.
"""
struct GrdDimMetadata{K,V} <: DimMetadata{K,V}
    val::Dict{K,V}
end
GrdDimMetadata() = GrdDimMetadata(Dict())


# Grd attributes wrapper

struct GrdAttrib{T,F,A}
    filename::F
    attrib::A
end
GrdAttrib(filename::AbstractString) = begin
    filename = first(splitext(filename))
    lines = readlines(filename * ".grd")
    entries = filter!(x -> !isempty(x) && !(x[1] == '['), lines)
    attrib = Dict(Pair(string.(strip.(match(r"([^=]+)=(.*)", st).captures[1:2]))...) for st in entries)
    T = datatype_translation[attrib["datatype"]]

    GrdAttrib{T,typeof(filename),typeof(attrib)}(filename, attrib)
end

filename(grd::GrdAttrib) = grd.filename

dims(grd::GrdAttrib, usercrs=nothing) = begin
    attrib = grd.attrib
    crs = ProjString(grd.attrib["projection"])

    ncols, nrows, nbands = size(grd)

    xbounds = parse.(Float64, (attrib["xmin"], attrib["xmax"]))
    ybounds = parse.(Float64, (attrib["ymin"], attrib["ymax"]))

    xspan = (xbounds[2] - xbounds[1]) / ncols
    yspan = (ybounds[2] - ybounds[1]) / nrows

    # Not fully implemented yet
    latlon_metadata = GrdDimMetadata(Dict())

    latmode = ProjectedIndex(
        order=Ordered(Forward(), Reverse(), Reverse()), 
        span=Regular(xspan), 
        sampling=Intervals(Start()), 
        crs=crs, 
        usercrs=usercrs,
    )
    lat = Lat(LinRange(ybounds[1], ybounds[2] - yspan, nrows), latmode, latlon_metadata)
    lonmode = ProjectedIndex(
        order=Ordered(),
        span=Regular(yspan), 
        sampling=Intervals(Start()), 
        crs=crs, 
        usercrs=usercrs,
    ) 
    lon = Lon(LinRange(xbounds[1], xbounds[2] - xspan, ncols), lonmode, latlon_metadata)
    band = Band(1:nbands; mode=Categorical(Ordered()))
    lon, lat, band
end

metadata(grd::GrdAttrib) = begin
    metadata = GrdMetadata() 
    for key in ("creator", "created", "history")
        val = get(grd.attrib, key, "")
        if val != ""
            metadata[key] = val
        end
    end
    metadata
end

missingval(grd::GrdAttrib{T}) where T = parse(T, grd.attrib["nodatavalue"])

name(grd::GrdAttrib) = get(grd.attrib, "layername", "")


Base.eltype(::GrdAttrib{T}) where T = T

Base.size(grd::GrdAttrib) = begin
    ncols = parse(Int, grd.attrib["ncols"])
    nrows = parse(Int, grd.attrib["nrows"])
    nbands = parse(Int, grd.attrib["nbands"])
    # The order is backwards to jula array order
    ncols, nrows, nbands
end

Base.Array(grd::GrdAttrib) = mmapgrd(Array, grd)


# Array ########################################################################

"""
    GrdArray(filename::String; refdims=(), name=nothing, window=(), usercrs=nothing)

An [`AbstractGeoArray`](@ref) that loads .grd files lazily from disk.

## Arguments
- `filename`: `String` pointing to a grd file. Extension is optional.

## Keyword arguments
- `name`: Name for the array. Will be loaded from `layername` if not supplied.
- `refdims`: Add dimension position array was sliced from. Mostly used programatically.
- `usercrs`: can be any CRS `GeoFormat` form GeoFormatTypes.jl, such as `WellKnownText`
- `window`: `Tuple` of `Dimension`, `Selector` or regular index to be applied when 
  loading the array. Can save on disk load time for large files.

# Example
```julia
array = GrdArray("folder/file.grd"; usercrs=EPSG(4326))
# Select Australia using 4326 coords, whatever the crs is underneath.
array[Lat(Between(-10, -43), Lon(113, 153))
```
"""
struct GrdArray{T,N,A,D<:Tuple,R<:Tuple,Na<:AbstractString,Me,Mi,W,S
               } <: DiskGeoArray{T,N,D,LazyArray{T,N}}
    filename::A
    dims::D
    refdims::R
    name::Na
    metadata::Me
    missingval::Mi
    window::W
    size::S
end
GrdArray(filename::String; kwargs...) = GrdArray(GrdAttrib(filename), kwargs...)
GrdArray(grd::GrdAttrib; refdims=(), name=nothing, window=(), usercrs=nothing) = begin
    attrib = grd.attrib

    dims_ = dims(grd, usercrs)
    metadata_ = metadata(grd)
    missingval_ = missingval(grd)
    size_ = map(length, dims_)
    if name isa Nothing
        name = GeoData.name(grd)
    end

    T = eltype(grd)
    N = length(size_)

    GrdArray{T,N,typeof.((grd.filename, dims_,refdims,name,metadata_,missingval_,window,size_))...
            }(grd.filename, dims_, refdims, name, metadata_, missingval_, window, size_)
end

# AbstractGeoArray methods

data(A::GrdArray) =
    mmapgrd(A) do mmap
        _window = maybewindow2indices(mmap, dims(A), window(A))
        readwindowed(mmap, _window)
    end

# Base methods

Base.getindex(A::GrdArray, I::Vararg{<:Union{<:Integer,<:AbstractArray}}) =
    mmapgrd(A) do mmap
        _window = maybewindow2indices(mmap, dims(A), window(A))
        _dims, _refdims = slicedims(slicedims(dims(A), refdims(A), _window)..., I)
        data = readwindowed(mmap, _window, I...)
        rebuild(A, data, _dims, _refdims)
    end
Base.getindex(A::GrdArray, i1::Integer, I::Vararg{<:Integer}) =
    mmapgrd(A) do mmap
        _window = maybewindow2indices(mmap, dims(A), window(A))
        readwindowed(mmap, _window, i1, I...)
    end

"""
    Base.write(filename::AbstractString, ::Type{GrdArray}, s::AbstractGeoArray)

Write a [`GrdArray`](@ref) to a .grd file, with a .gri header file. The extension of 
`filename` will be ignored. 

Currently the `metadata` field is lost on `write`. 
"""
Base.write(filename::String, ::Type{GrdArray}, A::AbstractGeoArray) = begin
    if hasdim(A, Band)
        correctedA = permutedims(A, (Lon, Lat, Band))# |>
            #a -> reorderindex(a, Forward()) |>
            #a -> reorderrelation(a, (Lon(Forward()), Lat(Reverse()), Band(Forward())))
        nbands = length(val(dims(A, Band)))
    else
        correctedA = permutedims(A, (Lon, Lat))# |>
            #a -> reorderindex(a, Forward()) |>
            #a -> reorderrelation(a, (Lon(Forward()), Lat(Reverse())))
        nbands = 1
    end
    # Remove extension
    filename = splitext(filename)[1]
    ncols, nrows = size(A)
    xmin, xmax = bounds(dims(A, Lon))
    ymin, ymax = bounds(dims(A, Lat))
    proj = convert(String, crs(dims(A, Lat)))
    datatype = rev_datatype_translation[eltype(A)]
    nodatavalue = missingval(A)
    minvalue = minimum(filter(x -> x != missingval(A), data(A)))
    maxvalue = maximum(filter(x -> x != missingval(A), data(A)))

    lat_array_ord = arrayorder(correctedA, Lat)
    if !(lat_array_ord isa Reverse)
        @warn "Array data order for Lat is `$lat_array_ord`, usualy `Reverse()`"
    end

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
    return
end


# AbstractGeoStack methods

GrdStack(filename; kwargs...) = DiskStack(filename; childtype=GrdArray, kwargs...)

querychild(f, ::Type{<:GrdArray}, filename::AbstractString, stack) = 
    f(GrdAttrib(filename))


# Utils ########################################################################
mmapgrd(f, A::GrdArray, mode...) = mmapgrd(f, filename(A), eltype(A), size(A), mode...) 
mmapgrd(f, grd::GrdAttrib, mode...) = mmapgrd(f, filename(grd), eltype(grd), size(grd), mode...) 
mmapgrd(f, filename::AbstractString, T::Type, size::Tuple, mode::AbstractString="r") = begin
    open(filename * ".gri", mode) do io
        mmap = Mmap.mmap(io, Array{T,length(size)}, size)
        output = f(mmap)
        close(io)
        output
    end
end

const datatype_translation = Dict{String, DataType}(
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

const rev_datatype_translation = Dict{DataType, String}(v => k for (k,v) in datatype_translation)
