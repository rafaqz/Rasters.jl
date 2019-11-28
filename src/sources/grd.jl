export GrdArray

struct GrdMetadata{M} <: AbstractArrayMetadata
    val::M 
end

struct GrdDimMetadata{M} <: AbstractDimMetadata
    val::M 
end

# Array ########################################################################

struct GrdArray{T,N,A,D<:Tuple,R<:Tuple,Me,Mi,Na,P,W,S} <: DiskGeoArray{T,N,D}
    filename::A
    dims::D
    refdims::R
    metadata::Me
    missingval::Mi
    name::Na
    window::W
    size::S
end

GrdArray(filepath::String; refdims=(), metadata=Dict(), window=()) = begin
    filepath = joinpath(dirname(filepath), first(splitext(filepath)))
    lines = readlines(filepath * ".grd")
    entries = filter!(x -> !isempty(x) && !(x[1] == '['), lines)
    data = Dict((c = match(r"([^=]+)=(.*)", st); string(c.captures[1]) => string(strip(c.captures[2]))) for st in entries)

    nrows = parse(Int, data["nrows"])
    ncols = parse(Int, data["ncols"])
    nbands = parse(Int, data["nbands"])

    T = datatype_translation[data["datatype"]]
    xmin, xmax = parse.(Float64, (data["xmin"], data["xmax"]))
    ymin, ymax = parse.(Float64, (data["ymin"], data["ymax"]))
    cellx = (xmax - xmin) / nrows
    celly = (ymax - ymin) / ncols

    lon = Lon(LinRange(ymin, ymax - celly, ncols); 
              grid=RegularGrid(; locus=Start(), span=celly))
    lat = Lat(LinRange(xmin, xmax - cellx, nrows); 
              grid=RegularGrid(; order=Ordered(Forward(), Reverse), locus=Start(), span=cellx)) 
    band = Band(1:nbands; grid=CategoricalGrid())
    dims = lon, lat, band
    for key in ("creator", "created", "history")
        val = get(data, key, nothing)
        if val !== nothing
            metadata[key] = history
        end
    end
    metadata = GrdMetadata(metadata)
    missingval = parse(T, data["nodatavalue"])
    name = get(data, "layername", "unnamed")

    _size = ncols, nrows, nbands

    GrdArray{T,length(_size),typeof.((filepath,dims,refdims,metadata,missingval,name))
            }(filepath, dims, refdims, metadata, missingval, name, window, _size)
end

Base.parent(A::GrdArray) =
    grdapply(mmap -> (w = windoworempty(A); w == () ? Array(mmap) : mmap[w...]), A)

Base.getindex(A::GrdArray, I::Vararg{<:Union{<:Integer,<:AbstractArray}}) = begin
    I = applywindow(A, I)
    rebuildsliced(A, grdapply(mmap -> mmap[I...], A), I)
end
Base.getindex(A::GrdArray, I::Vararg{<:Integer}) = 
    grdapply(mmap -> mmap[applywindow(A, I)...], A)

Base.setindex!(A::GrdArray, x, I::Vararg{<:Union{<:Integer,<:AbstractArray}}) = begin
    I = applywindow(A, I)
    grdapply(A -> A[I...] = x, A, "w+")
end

grdapply(f, A, mode="r") = begin
    open(filename(A) * ".gri", mode) do io
        mmap = Mmap.mmap(io, Array{eltype(A),3}, size(A))
        output = f(mmap)
        close(io)
        output 
    end
end

Base.write(filename::String, GrdArray, A::AbstractGeoArray) = begin
    # grid(dims(A) <: RegularGrid || throw(ArgumentError("Can only save `RegularGrid` arrays to a grd file"))
    # Remove extension
    filename = splitext(filename)[1]
    # Standardise dimensions
    # A = permutedims(A, [Lon, Lat])
    ncols, nrows = size(A)
    xmin, xmax = bounds(dims(A, Lon()))
    ymin, ymax = bounds(dims(A, Lat()))
    proj = convert(String, projection(A))
    datatype = rev_datatype_translation[eltype(A)]
    nodatavalue = missingval(A)
    minvalue = minimum(filter(x -> x != missingval(A), parent(A)))
    maxvalue = maximum(filter(x -> x != missingval(A), parent(A)))
    nbands = hasdim(A, Band()) ? length(val(dims(A, Band()))) : 1

    # Data: gri file
    open(filename  * ".gri", "w") do IO
        write(IO, parent(A))
    end

    # Metadata: grd file
    open(filename * ".grd", "w") do IO
        write(IO,
            """
            [general]
            creator=GeoData.jl
            created= $(string(now()))
            [georeference]
            nrows= $(nrows)
            ncols= $(ncols)
            xmin= $(xmin)
            ymin= $(ymin)
            xmax= $(xmax)
            ymax= $(ymax)
            projection= $(proj)
            [data]
            datatype= $(datatype)
            nodatavalue= $(nodatavalue)
            byteorder= little
            nbands= $nbands
            minvalue= $(minvalue)
            maxvalue= $(maxvalue)
            [description]
            layername= $(name(A))
            """
        )
    end
    return
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
