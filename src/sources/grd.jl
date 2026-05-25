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
    attrib = Dict(pairs)
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
    crs = if crs isa NoKW 
        str = attrib["projection"]
        str == "" ? nothing : ProjString(str)
    else
        crs
    end
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

function missingval(A::RasterDiskArray{GRDsource,T}, args...) where T
    _grd_mv(T, A.attrib.attrib)
end

function _grd_mv(::Type{T}, md; verbose=true) where T
    if haskey(md, "nodatavalue")
        ndv = md["nodatavalue"]
        ndv === "nothing" && return nothing
        try
            return parse(T, ndv)
        catch
            verbose && @warn "nodatavalue $(ndv) is not convertible to data type $T. `missingval` set to `nothing`."
            return nothing
        end
    else
        return nothing
    end
end

_sizeof(A::GRDdataset{T}) where T = sizeof(T) * prod(size(A))
_sizeof(A::RasterDiskArray{GRDsource}) = _sizeof(A.attrib)

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
    verbose=true,
    write=true,
    missingval=nokw,
    chunks=nokw,
    scale=nokw,
    offset=nokw,
    coerce=nokw,
    eltype=Missings.nonmissingtype(eltype(A)),
    f=identity,
    kw...
)
    check_can_write(filename, force)
    write = f === identity ? write : true
    haskey(REVGRD_DATATYPE_TRANSLATION, eltype) || throw(ArgumentError("""
       Element type $eltype cannot be written to grd file. Convert it to one of $(keys(REVGRD_DATATYPE_TRANSLATION)),
       usually by broadcasting the desired type constructor over the `Raster`, e.g. `newrast = Float32.(rast)`"))
       """
    ))
    isnokwornothing(scale) && isnokwornothing(offset) || throw(ArgumentError("Cant write scale or offset to .grd files"))
    chunks isa NoKW || @warn "specifying chunks not supported for .grd files"
    # Missing values
    missingval_pair = _write_missingval_pair(A, missingval; eltype, verbose)

    # Missing values

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

    # Data: write a raw gri file from the array
    mod = _mod(eltype, missingval_pair, scale, offset, coerce)
    gri_filename = filename * ".gri"
    _write_gri(gri_filename, Val{source_eltype(mod)}(), mod, parent(correctedA))
    _write_grd(filename, eltype, dims(A), missingval_pair[1], name(A))

    if write
        _mmapgrd(filename, source_eltype(mod), size(A); write=true) do M
            f(_maybe_modify(M, mod))
        end
    end

    return filename * ".grd"
end

function _write_gri(filename, v, ::NoMod, A::Array{T}) where T
    open(filename; write=true, lock=false) do io
        write(io, A)
    end
end
function _write_gri(filename, v, mod, A::AbstractArray)
    open(filename; write=true, lock=false) do io
        # Avoid `Ref` allocations
        ref = Ref{source_eltype(mod)}(_invertmod(v, first(A), mod))
        for x in A # We are modifying the source array so invert the modifications
            ref[] = _invertmod(v, x, mod)
            write(io, ref)
        end
    end
end

function _write_grd(filename, T, dims, missingval, name)
    filename = splitext(filename)[1]

    x, y = map(DD.dims(dims, (X(), Y()))) do d
        convertlookup(Projected, d)
    end

    nbands = hasdim(dims, Band) ? length(DD.dims(dims, Band)) : 1
    ncols, nrows = length(x), length(y)
    xmin, xmax = bounds(x)
    ymin, ymax = bounds(y)
    proj = isnothing(crs(x)) ? "" : convert(String, convert(ProjString, crs(x)))
    datatype = REVGRD_DATATYPE_TRANSLATION[T]
    nodatavalue = missingval

    # Metadata: grd file
    open(filename * ".grd"; write=true, lock=false) do IO
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
            [description]
            layername= $name
            """
        )
    end
end

# Rasters methods
function _open(f, source::GRDsource, filename::AbstractString; kw...)
    isfile(filename) || _filenotfound_error(filename)
    _open(f, source, GRDdataset(filename); kw...)
end
function _open(f, source::GRDsource, ds::GRDdataset; write=false, kw...)
    _mmapgrd(ds; write) do mm
        A = RasterDiskArray{GRDsource}(mm, DA.eachchunk(ds), DA.haschunks(ds), ds)
        _open(f, source, A; kw...)
    end
end
_open(f, ::GRDsource, A::RasterDiskArray; mod=NoMod(), kw...) =
    cleanreturn(f(_maybe_modify(A, mod)))

haslayers(::GRDsource) = false

Raster(ds::RasterDiskArray{<:GRDsource}; kw...) = _raster(ds; kw...)


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
