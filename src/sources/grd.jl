const GRD_X_ORDER = ForwardOrdered()
const GRD_Y_ORDER = ReverseOrdered()
const GRD_BAND_ORDER = ForwardOrdered()

const GRD_DATATYPE_TRANSLATION = Dict{String, DataType}(
    "LOG1S" => Bool,
    "INT1S" => Int8,
    "INT2S" => Int16,
    "INT4S" => Int32,
    "INT8S" => Int64,
    "INT1U" => UInt8,
    "INT2U" => UInt16,
    "INT4U" => UInt32,
    "INT8U" => UInt64,
    "FLT4S" => Float32,
    "FLT8S" => Float64
)
const REVGRD_DATATYPE_TRANSLATION =
    Dict{DataType, String}(v => k for (k,v) in GRD_DATATYPE_TRANSLATION)

# GRD attributes wrapper. Only used during file load, for dispatch.
struct GRDdataset{T,F}
    filename::F
    attrib::Dict{String,String}
    write::Bool
end
function GRDdataset(filename::AbstractString; write=false)
    filename = first(splitext(filename))
    lines = readlines(filename * ".grd")
    entries = filter!(x -> !isempty(x) && !(x[1] == '['), lines)
    matches = (match(r"([^=]+)=(.*)", st) for st in entries)
    captures = (string.(strip.(m.captures[1:2])) for m in matches)
    pairs = map(c -> Pair(c[1], c[2]), captures)
    attrib = Dict(pairs...)
    T = GRD_DATATYPE_TRANSLATION[attrib["datatype"]]
    GRDdataset{T,typeof(filename)}(filename, attrib, write)
end

attrib(grd::GRDdataset) = grd.attrib
filename(grd::GRDdataset) = grd.filename
filekey(grd::GRDdataset, name::NoKW) = get(attrib(grd), "layername", Symbol(""))
filekey(A::RasterDiskArray{GRDsource}, name) = filekey(A.attrib, name)

Base.eltype(::GRDdataset{T}) where T = T
function Base.size(grd::GRDdataset)
    ncols = parse(Int, grd.attrib["ncols"])
    nrows = parse(Int, grd.attrib["nrows"])
    nbands = parse(Int, grd.attrib["nbands"])
    # The order is backwards to jula array order
    ncols, nrows, nbands
end
Base.Array(grd::GRDdataset) = _mmapgrd(Array, grd)

function DiskArrays.eachchunk(A::GRDdataset) 
    size_ = size(A)
    DiskArrays.GridChunks(size_, size_)
end
DiskArrays.haschunks(::GRDdataset) = DiskArrays.Unchunked()

function _dims(A::RasterDiskArray{GRDsource}, crs=nokw, mappedcrs=nokw)
    attrib = A.attrib.attrib
    crs = crs isa NoKW ? ProjString(attrib["projection"]) : crs
    mappedcrs = mappedcrs isa NoKW ? nothing : mappedcrs

    xsize, ysize, nbands = size(A)

    xbounds = parse.(Float64, (attrib["xmin"], attrib["xmax"]))
    ybounds = parse.(Float64, (attrib["ymin"], attrib["ymax"]))

    # Always intervals
    xspan = (xbounds[2] - xbounds[1]) / xsize
    yspan = (ybounds[1] - ybounds[2]) / ysize

    # Not fully implemented yet
    xy_metadata = _metadatadict(GRDsource())

    xindex = range(; start=xbounds[1], stop=xbounds[2] - xspan, length=xsize)
    yindex = range(; start=ybounds[2] + yspan, stop=ybounds[1], length=ysize)

    xlookup = Projected(xindex;
        order=GRD_X_ORDER,
        span=Regular(xspan),
        sampling=Intervals(Start()),
        metadata=xy_metadata,
        crs=crs,
        mappedcrs=mappedcrs,
        dim=X()
    )
    ylookup = Projected(yindex;
        order=GRD_Y_ORDER,
        span=Regular(yspan),
        sampling=Intervals(Start()),
        metadata=xy_metadata,
        crs=crs,
        mappedcrs=mappedcrs,
        dim=Y()
    )

    x = X(xlookup)
    y = Y(ylookup)
    band = Band(Categorical(1:nbands; order=GRD_BAND_ORDER))
    return x, y, band
end

_name(A::RasterDiskArray{GRDsource}) = Symbol(get(A.attrib.attrib, "layername", ""))

function _metadata(A::RasterDiskArray{GRDsource}, args...)
    metadata = _metadatadict(GRDsource())
    for key in ("creator", "created", "history")
        val = get(A.attrib.attrib, key, "")
        if val != ""
            metadata[key] = val
        end
    end
    metadata
end

function missingval(A::RasterDiskArray{GRDsource,T}) where T
    if haskey(A.attrib.attrib, "nodatavalue")
        ndv = A.attrib.attrib["nodatavalue"]
        ndv === "nothing" && return nothing
        try
            return parse(T, ndv)
        catch
            @warn "nodatavalue $(ndv) is not convertible to data type $T. `missingval` set to `nothing`."
            return nothing
        end
    else
        return nothing
    end
end

_sizeof(A::GRDdataset{T}) where T = sizeof(T) * prod(size(A))
_sizeof(A::RasterDiskArray{GRDsource}) = _sizeof(A.attrib)


# Array ########################################################################

function FileArray{GRDsource}(A::RasterDiskArray{<:Any,T}, filename=filename(A.attrib); kw...) where T
    filename = first(splitext(filename))
    eachchunk = DiskArrays.eachchunk(A)
    haschunks = DiskArrays.haschunks(A)
    FileArray{GRDsource,T,3}(filename, size(A); eachchunk, haschunks, kw...)
end

# Base methods

"""
    Base.write(filename::AbstractString, ::Type{GRDsource}, s::AbstractRaster; kw...)

Write a `Raster` to a .grd file with a .gri header file.

This method is called automatically if you `write` a `Raster` 
with a `.grd` or `.gri` extension. 

## Keywords 

$FORCE_KEYWORD

If this method is called directly the extension of `filename` will be ignored.

Returns the base of `filename` with a `.grd` extension.
"""
function Base.write(filename::String, ::GRDsource, A::AbstractRaster; 
    force=false, 
    missingval=nokw,
    chunks=nokw,
    kw...
)
    chunks isa NoKW || @warn "specifying chunks not supported for .grd files"
    check_can_write(filename, force)
    A = _maybe_use_type_missingval(A, GRDsource(), missingval)
    if hasdim(A, Band)
        correctedA = permutedims(A, (X, Y, Band)) |>
            a -> reorder(a, (X(GRD_X_ORDER), Y(GRD_Y_ORDER), Band(GRD_BAND_ORDER)))
        nbands = length(val(dims(correctedA, Band)))
    else
        correctedA = permutedims(A, (X, Y)) |>
            a -> reorder(a, (X(GRD_X_ORDER), Y(GRD_Y_ORDER)))
        nbands = 1
    end
    # Remove extension
    filename = splitext(filename)[1]
    minvalue = minimum(skipmissing(A))
    maxvalue = maximum(skipmissing(A))
    _write_grd(filename, eltype(A), dims(A), Rasters.missingval(A), minvalue, maxvalue, name(A))

    # Data: gri file
    open(filename * ".gri", write=true) do IO
        write(IO, parent(correctedA))
    end

    return filename * ".grd"
end

function _write_grd(filename, T, dims, missingval, minvalue, maxvalue, name)
    filename = splitext(filename)[1]

    x, y = map(DD.dims(dims, (X(), Y()))) do d
        convertlookup(Projected, d)
    end

    nbands = hasdim(dims, Band) ? length(DD.dims(dims, Band)) : 1
    ncols, nrows = length(x), length(y)
    xmin, xmax = bounds(x)
    ymin, ymax = bounds(y)
    proj = convert(String, convert(ProjString, crs(x)))
    datatype = REVGRD_DATATYPE_TRANSLATION[T]
    nodatavalue = missingval

    # Metadata: grd file
    open(filename * ".grd"; write=true) do IO
        write(IO,
            """
            [general]
            creator=Rasters.jl
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
            layername= $name
            """
        )
    end
end

function create(filename, ::GRDsource, T::Type, dims::DD.DimTuple; 
    name="layer", metadata=nothing, missingval=nothing, lazy=true, 
)
    # Remove extension
    basename = splitext(filename)[1]
    minvalue = maxvalue = zero(T)
    sze = map(length, DD.dims(dims, (XDim, YDim, Band)))

    # Metadata: grd file
    _write_grd(basename, T, dims, missingval, minvalue, maxvalue, name)

    # Data: gri file
    open(basename * ".gri", write=true) do IO
        write(IO, FillArrays.Zeros(sze))
    end
    return Raster(filename; source=GRDsource(), lazy)
end

# AbstractRasterStack methods

# Custom `open` because the data and metadata objects are separate
# Here we _mmapgrd instead of `_open`
function Base.open(f::Function, A::FileArray{GRDsource}, args...; write=A.write)
    _mmapgrd(mm -> f(RasterDiskArray{GRDsource}(mm, A.eachchunk, A.haschunks)), A; write)
end

function _open(f, ::GRDsource, filename::AbstractString; write=false, name=nokw, group=nokw)
    isfile(filename) || _filenotfound_error(filename)
    attr = GRDdataset(filename)
    _mmapgrd(attr; write) do mm
        A = RasterDiskArray{GRDsource}(mm, DA.eachchunk(attr), DA.haschunks(attr), attr)
        f(A)
    end
end
_open(f, ::GRDsource, attrib::GRDdataset; kw...) = f(attrib)

haslayers(::GRDsource) = false

# Utils ########################################################################

function _mmapgrd(f, x::Union{FileArray,GRDdataset}; kw...)
    _mmapgrd(f, filename(x), eltype(x), size(x); kw...)
end
function _mmapgrd(f, filename::AbstractString, T::Type, size::Tuple; write=false)
    arg = write ? "r+" : "r"  
    basename = splitext(filename)[1]
    open(basename * ".gri", arg) do io
        mmap = Mmap.mmap(io, Array{T,length(size)}, size)
        output = f(mmap)
        close(io)
        output
    end
end

# precompilation
# function _precompile(::Type{GRDsource})
#     ccall(:jl_generating_output, Cint, ()) == 1 || return nothing

#     T = UInt16
#     for T in (Any, UInt8, UInt16, Int16, UInt32, Int32, Int64, Float32, Float64)
#         precompile(GRDdataset, (String,))
#         DS = Rasters.GRDdataset{T,String}
#         precompile(crs, (DS,))
#         precompile(Rasters.FileArray, (DS, String))
#         precompile(dims, (DS,))
#         precompile(dims, (DS,WellKnownText{GeoFormatTypes.CRS},Nothing))
#         precompile(dims, (DS,WellKnownText{GeoFormatTypes.CRS},EPSG))
#         precompile(dims, (DS,WellKnownText{GeoFormatTypes.CRS},ProjString))
#         precompile(dims, (DS,WellKnownText{GeoFormatTypes.CRS},WellKnownText{GeoFormatTypes.CRS}))
#         precompile(metadata, (DS, ))
#         precompile(metadata, (DS, Symbol))
#         precompile(missingval, (DS,))
#         precompile(Raster, (DS, String, Nothing))
#         precompile(Raster, (DS, String, Symbol))
#     end
# end

# _precompile(GRDsource)
