export GrdArray, GrdMetadata, GrdDimMetadata


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
GrdArray(filename::String; refdims=(), name=nothing, 
         metadata=GrdMetadata(), window=(), usercrs=nothing) = begin
    filename = first(splitext(filename))
    lines = readlines(filename * ".grd")
    entries = filter!(x -> !isempty(x) && !(x[1] == '['), lines)
    data = Dict((c = match(r"([^=]+)=(.*)", st); string(c.captures[1]) => string(strip(c.captures[2]))) for st in entries)

    nrows = parse(Int, data["nrows"])
    ncols = parse(Int, data["ncols"])
    nbands = parse(Int, data["nbands"])
    _size = ncols, nrows, nbands

    T = datatype_translation[data["datatype"]]
    N = length(_size)
    xmin, xmax = parse.(Float64, (data["xmin"], data["xmax"]))
    ymin, ymax = parse.(Float64, (data["ymin"], data["ymax"]))
    crs = ProjString(data["projection"])
    cellx = (xmax - xmin) / nrows
    celly = (ymax - ymin) / ncols
    # Not fully implemented yet
    latlon_metadata = GrdDimMetadata(Dict())

    latmode = ProjectedIndex(
        order=Ordered(Forward(), Reverse(), Forward()), 
        span=Regular(cellx), 
        sampling=Intervals(Start()), 
        crs=crs, 
        usercrs=usercrs,
    )
    lat = Lat(LinRange(xmin, xmax - cellx, nrows), latmode, latlon_metadata)
    lonmode = ProjectedIndex(
        order=Ordered(),
        span=Regular(celly), 
        sampling=Intervals(Start()), 
        crs=crs, 
        usercrs=usercrs,
    ) 
    lon = Lon(LinRange(ymin, ymax - celly, ncols), lonmode, latlon_metadata)
    band = Band(1:nbands; mode=Categorical(Ordered()))
    dims = lon, lat, band
    for key in ("creator", "created", "history")
        val = get(data, key, "")
        if val != ""
            metadata[key] = val
        end
    end
    missingval = parse(T, data["nodatavalue"])
    if name isa Nothing
        get(data, "layername", "")
    end

    GrdArray{T,N,typeof.((filename,dims,refdims,name,metadata,missingval,window,_size))...
            }(filename, dims, refdims, name, metadata, missingval, window, _size)
end

# AbstractGeoStack methods

data(A::GrdArray) =
    grdapply(A) do mmap
        _window = maybewindow2indices(mmap, dims(A), window(A))
        readwindowed(mmap, _window)
    end

# Base methods

Base.getindex(A::GrdArray, I::Vararg{<:Union{<:Integer,<:AbstractArray}}) =
    grdapply(A) do mmap
        _window = maybewindow2indices(mmap, dims(A), window(A))
        _dims, _refdims = slicedims(slicedims(dims(A), refdims(A), _window)..., I)
        data = readwindowed(mmap, _window, I...)
        rebuild(A, data, _dims, _refdims)
    end
Base.getindex(A::GrdArray, i1::Integer, I::Vararg{<:Integer}) =
    grdapply(A) do mmap
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
        correctedA = permutedims(A, (Lon, Lat, Band)) |>
            a -> reorderindex(a, Forward()) |>
            a -> reorderrelation(a, Forward())
        nbands = length(val(dims(A, Band)))
    else
        correctedA = permutedims(A, (Lon, Lat)) |>
            a -> reorderindex(a, Forward()) |>
            a -> reorderrelation(a, Forward())
        nbands = 1
    end
    # Remove extension
    filename = splitext(filename)[1]
    ncols, nrows = size(A)
    xmin, xmax = bounds(dims(A, Lat))
    ymin, ymax = bounds(dims(A, Lon))
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
    open(filename  * ".gri", "w") do IO
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


# Utils ########################################################################

grdapply(f, A, mode="r") = begin
    open(filename(A) * ".gri", mode) do io
        mmap = Mmap.mmap(io, Array{eltype(A),3}, size(A))
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
