export GrdArray, GrdMetadata, GrdDimMetadata


# Metadata ########################################################################

struct GrdMetadata{K,V} <: ArrayMetadata{K,V}
    val::Dict{K,V}
end

struct GrdDimMetadata{K,V} <: DimMetadata{K,V}
    val::Dict{K,V}
end


# Array ########################################################################

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

GrdArray(filepath::String; 
         refdims=(), 
         name=nothing, 
         metadata=GrdMetadata(Dict()), 
         window=(),
         selectorcrs=nothing,
        ) = begin
    filepath = first(splitext(filepath))
    lines = readlines(filepath * ".grd")
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
    latlon_metadata = GrdDimMetadata(Dict(:crs => crs))

    lat = Lat(LinRange(xmin, xmax - cellx, nrows);
              grid=ProjectedGrid(; order=Ordered(Forward(), Reverse(), Forward()), 
                               locus=Start(), step=cellx, crs=crs, selectorcrs=selectorcrs),
              metadata=latlon_metadata)
    lon = Lon(LinRange(ymin, ymax - celly, ncols);
              grid=ProjectedGrid(; locus=Start(), step=celly, crs=crs, selectorcrs=selectorcrs),
              metadata=latlon_metadata)
    band = Band(1:nbands; grid=CategoricalGrid())
    dims = lon, lat, band
    for key in ("creator", "created", "history")
        val = get(data, key, "")
        if val != ""
            metadata[key] = val
        end
    end
    missingval = parse(T, data["nodatavalue"])
    if !(name isa String)
        name = get(data, "layername", "")
    end

    GrdArray{T,N,typeof.((filepath,dims,refdims,name,metadata,missingval,window,_size))...
            }(filepath, dims, refdims, name, metadata, missingval, window, _size)
end

data(A::GrdArray) =
    grdapply(A) do mmap
        _window = maybewindow2indices(mmap, dims(A), window(A))
        readwindowed(mmap, _window)
    end
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


# Base.setindex!(A::GrdArray, x, I::Vararg{<:Union{<:Integer,<:AbstractArray}}) =
    # grdapply(mmap -> mmap[applywindow(A, I)...] = x, A, "w+")

Base.write(filename::String, ::Type{GrdArray}, A::AbstractGeoArray) = begin
    # grid(dims(A) <: RegularGrid || throw(ArgumentError("Can only save `RegularGrid` arrays to a grd file"))
    # Remove extension
    filename = splitext(filename)[1]
    # Standardise dimensions
    # A = permutedims(A, [Lon, Lat])
    ncols, nrows = size(A)
    xmin, xmax = bounds(dims(A, Lat()))
    ymin, ymax = bounds(dims(A, Lon()))
    proj = convert(String, crs(dims(A, Lat())))
    datatype = rev_datatype_translation[eltype(A)]
    nodatavalue = missingval(A)
    minvalue = minimum(filter(x -> x != missingval(A), data(A)))
    maxvalue = maximum(filter(x -> x != missingval(A), data(A)))
    nbands = hasdim(A, Band()) ? length(val(dims(A, Band()))) : 1

    # Data: gri file
    open(filename  * ".gri", "w") do IO
        write(IO, data(A))
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
