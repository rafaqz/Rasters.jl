const AG = ArchGDAL

const GDAL_X_ORDER = ForwardOrdered()
const GDAL_Y_ORDER = ReverseOrdered()

const GDAL_X_LOCUS = Start()
const GDAL_Y_LOCUS = Start()

# Array ########################################################################

function FileArray(raster::AG.RasterDataset{T}, filename; kw...) where {T}
    eachchunk, haschunks = DA.eachchunk(raster), DA.haschunks(raster)
    FileArray{GDALfile,T,3}(filename, size(raster); eachchunk, haschunks, kw...)
end

cleanreturn(A::AG.RasterDataset) = Array(A)

haslayers(::Type{GDALfile}) = false

"""
    Base.write(filename::AbstractString, ::Type{GDALfile}, A::AbstractRaster; kw...)

Write a `Raster` to file using GDAL.

# Keywords

- `driver::String`: a GDAL driver name. Guessed from the filename extension by default.
- `compress::String`: GeoTIFF compression flag. "DEFLATE" by default.
- `tiled::Bool`: GeoTiff tiling. Defaults to `true`.

Returns `filename`.
"""
function Base.write(
    filename::AbstractString, ::Type{GDALfile}, A::AbstractRaster{T,2}; kw...
) where T
    all(hasdim(A, (X, Y))) || error("Array must have Y and X dims")

    correctedA = if lookup(A, X) isa AffineProjected
        A
    else
        _maybe_permute_to_gdal(A) |>
            a -> noindex_to_sampled(a) |>
            a -> reorder(a, (X(GDAL_X_ORDER), Y(GDAL_Y_ORDER)))
    end
    nbands = 1
    _gdalwrite(filename, correctedA, nbands; kw...)
end
function Base.write(
    filename::AbstractString, ::Type{GDALfile}, A::AbstractRaster{T,3}; kw...
) where T
    all(hasdim(A, (X, Y))) || error("Array must have Y and X dims")
    hasdim(A, Band()) || error("Must have a `Band` dimension to write a 3-dimensional array")

    correctedA = if lookup(A, X) isa AffineProjected
        A
    else
        _maybe_permute_to_gdal(A) |>
            a -> noindex_to_sampled(a) |>
            a -> reorder(a, (X(GDAL_X_ORDER), Y(GDAL_Y_ORDER)))
    end

    nbands = size(correctedA, Band())
    _gdalwrite(filename, correctedA, nbands; kw...)
end

function create(filename, ::Type{GDALfile}, T::Type, dims::DD.DimTuple;
    missingval=nothing, metadata=nothing, name=nothing, keys=(name,),
    driver=AG.extensiondriver(filename), compress="DEFLATE", chunk=nothing,
    parent=nothing, lazy=true, 
)
    if !(keys isa Nothing || keys isa Symbol) && length(keys) > 1
        throw(ArgumentError("GDAL cant write more than one layer per file, but keys $keys have $(length(keys))"))
    end
    x, y = map(DD.dims(dims, (XDim, YDim))) do d
        lookup(d) isa NoLookup ? set(d, Sampled) : d
    end
    x = reorder(x, GDAL_X_ORDER)
    y = reorder(y, GDAL_Y_ORDER)
    T = Missings.nonmissingtype(T)

    if ismissing(missingval)
        missingval = _writeable_missing(T)
    end
    if T === Int64
        @info "GDAL cannot create `Int64` files, using `Int32` instead"
        T = Int32
    end
    if missingval isa Int64
        missingval = Int32(missingval)
    end

    if hasdim(dims, Band)
        b = DD.dims(dims, Band)
        nbands = length(b)
        newdims = (x, y, b)
    else
        nbands = 1
        newdims = (x, y)
    end
    kw = (width=length(x), height=length(y), nbands=nbands, dtype=T)
    gdaldriver = AG.getdriver(driver)
    if driver == "GTiff"
        # TODO implement chunking
        options = ["COMPRESS=$compress", "TILED=NO"]
        AG.create(filename; driver=gdaldriver, options=options, kw...) do ds
            _gdalsetproperties!(ds, newdims, missingval)
        end
    else
        # Create a tif and copy it to `filename`, as ArchGDAL.create
        # does not support direct creation of ASCII etc. rasters
        ArchGDAL.create(tempname() * ".tif"; driver=AG.getdriver("GTiff"), kw...) do ds
            _gdalsetproperties!(ds, newdims, missingval)
            AG.copy(ds; filename=filename, driver=gdaldriver) |> AG.destroy
        end
    end
    if hasdim(dims, Band)
        return Raster(filename; source=GDALfile, lazy)
    else
        return view(Raster(filename; source=GDALfile, lazy), Band(1))
    end
end

# DimensionalData methods for ArchGDAL types ###############################

# We allow passing in crs and mappedcrs manually
function DD.dims(raster::AG.RasterDataset, crs=nothing, mappedcrs=nothing)
    gt = try
        AG.getgeotransform(raster)
    catch
        GDAL_EMPTY_TRANSFORM
    end
    xsize, ysize = size(raster)
    nbands = AG.nraster(raster)
    bandnames = _gdal_bandnames(raster, nbands)
    band = if all(==(""), bandnames)
        Band(Categorical(1:nbands; order=ForwardOrdered()))
    else
        Band(Categorical(bandnames; order=Unordered()))
    end

    crs = crs isa Nothing ? Rasters.crs(raster) : crs
    xy_metadata = metadata(raster)

    # Output Sampled index dims when the transformation is lat/lon alligned,
    # otherwise use Transformed index, with an affine map.
    if _isalligned(gt)
        xstep = gt[GDAL_WE_RES]
        if xstep > 0 
            xmin = gt[GDAL_TOPLEFT_X]
            xmax = gt[GDAL_TOPLEFT_X] + xstep * (xsize - 1)
            xorder = ForwardOrdered()
        else
            xmin = gt[GDAL_TOPLEFT_X] + xstep 
            xmax = gt[GDAL_TOPLEFT_X] + xstep * xsize
            xorder = ReverseOrdered()
        end
        xindex = LinRange(xmin, xmax, xsize)
        xorder = xstep > 0 ? ForwardOrdered() : ReverseOrdered()

        ystep = gt[GDAL_NS_RES] # A negative number
        if ystep > 0 
            ymax = gt[GDAL_TOPLEFT_Y]
            ymin = gt[GDAL_TOPLEFT_Y] + ystep * (ysize - 1)
            yorder = ForwardOrdered()
        else
            ymax = gt[GDAL_TOPLEFT_Y] + ystep
            ymin = gt[GDAL_TOPLEFT_Y] + ystep * ysize
            yorder = ReverseOrdered()
        end
        yindex = LinRange(ymax, ymin, ysize)

        # Spatial data defaults to area/inteval
        xsampling, ysampling = if _gdalmetadata(raster.ds, "AREA_OR_POINT") == "Point"
            Points(), Points()
        else
            # GeoTiff uses the "pixelCorner" convention
            Intervals(GDAL_X_LOCUS), Intervals(GDAL_Y_LOCUS)
        end

        xlookup = Projected(xindex;
            order=xorder,
            span=Regular(step(xindex)),
            sampling=xsampling,
            metadata=xy_metadata,
            crs=crs,
            mappedcrs=mappedcrs,
        )
        ylookup = Projected(yindex;
            order=yorder,
            sampling=ysampling,
            # Use the range step as is will be different to ystep due to float error
            span=Regular(step(yindex)),
            metadata=xy_metadata,
            crs=crs,
            mappedcrs=mappedcrs,
        )
        x = X(xlookup)
        y = Y(ylookup)

        DD.format((x, y, band), map(Base.OneTo, (xsize, ysize, nbands)))
    else
        affinemap = _geotransform2affine(gt)
        x = X(AffineProjected(affinemap; crs, mappedcrs, metadata=xy_metadata, dim=X()))
        y = Y(AffineProjected(affinemap; crs, mappedcrs, metadata=xy_metadata, dim=Y()))

        DD.format((x, y, band), map(Base.OneTo, (xsize, ysize, nbands)))
    end
end

function DD.metadata(raster::AG.RasterDataset, args...)
    band = AG.getband(raster.ds, 1)
    # color = AG.getname(AG.getcolorinterp(band))
    scale = AG.getscale(band)
    offset = AG.getoffset(band)
    # norvw = AG.noverview(band)
    path = first(AG.filelist(raster))
    units = AG.getunittype(band)
    upair = units == "" ? () : ("units"=>units,)
    _metadatadict(GDALfile, "filepath"=>path, "scale"=>scale, "offset"=>offset, upair...)
end

# Rasters methods for ArchGDAL types ##############################

# Create a Raster from a memory-backed dataset
Raster(ds::AG.Dataset; kw...) = Raster(AG.RasterDataset(ds); kw...)
function Raster(ds::AG.RasterDataset;
    crs=crs(ds), mappedcrs=nothing,
    dims=dims(ds, crs, mappedcrs),
    refdims=(), name=Symbol(""),
    metadata=metadata(ds),
    missingval=missingval(ds)
)
    args = dims, refdims, name, metadata, missingval
    filelist = AG.filelist(ds)
    if length(filelist) > 0
        filename = first(filelist)
        return Raster(FileArray(ds, filename), args...)
    else
        return Raster(Array(ds), args...)
    end
end

function missingval(raster::AG.RasterDataset, args...)
    # We can only handle data where all bands have the same missingval
    band = AG.getband(raster.ds, 1)
    nodata = AG.getnodatavalue(band)
    if nodata isa Nothing
        return nothing
    else
        return _gdalconvertmissing(eltype(band), nodata)
    end
end

function crs(raster::AG.RasterDataset, args...)
    WellKnownText(GeoFormatTypes.CRS(), string(AG.getproj(raster.ds)))
end

# ArchGDAL methods for Rasters types ####################################

# Extend ArchGDAL RasterDataset and Dataset to accept AbstractRaster as input
function AG.Dataset(f::Function, A::AbstractRaster; kw...)
    AG.RasterDataset(A; kw...) do rds
        f(rds.ds)
    end
end
function AG.RasterDataset(f::Function, A::AbstractRaster; filename=nothing)
    all(hasdim(A, (XDim, YDim))) || throw(ArgumentError("`AbstractRaster` must have both an `XDim` and `YDim` to use be converted to an ArchGDAL `Dataset`"))
    if ndims(A) === 3
        thirddim = otherdims(A, (X, Y))[1]
        thirddim isa Band || throw(ArgumentError("ArchGDAL can't handle $(basetypeof(thirddim)) dims - only XDim, YDim, and Band"))
    elseif ndims(A) > 3
        throw(ArgumentError("ArchGDAL can only accept 2 or 3 dimensional arrays"))
    end

    # block_x, block_y = DA.max_chunksize(DA.eachchunk(A))
    A_p = _maybe_permute_to_gdal(A)
    dataset = _unsafe_gdal_ds(A_p; filename)
    try
        rds = AG.RasterDataset(dataset)
        open(A_p) do a
            rds .= parent(a)
        end
        f(rds)
    finally
        AG.destroy(dataset)
    end
end


# Utils ########################################################################

function _open(f, ::Type{GDALfile}, filename::AbstractString; write=false, kw...)
    # Handle url filenames
    # /vsicurl/ is added to urls for GDAL, /vsimem/ for in memory
    if length(filename) >= 8 
        if (filename[1:7] == "http://" || filename[1:8] == "https://")
           filename = "/vsicurl/" * filename
        elseif !(filename[1:8] in ("/vsicurl", "/vsimem/"))
            # check the file actually exists because GDALs error is unhelpful
            isfile(filename) || _filenotfound_error(filename)
        end
    else
        isfile(filename) || _filenotfound_error(filename)
    end
    flags = write ? (; flags=AG.OF_UPDATE) : ()
    AG.readraster(cleanreturn ∘ f, filename; flags...)
end
_open(f, ::Type{GDALfile}, ds::AG.RasterDataset; kw...) = cleanreturn(f(ds))


_gdalconvertmissing(T::Type{<:AbstractFloat}, x::Real) = convert(T, x)
function _gdalconvertmissing(T::Type{<:Integer}, x::AbstractFloat)
    if trunc(x) === x
        convert(T, x)
    else
        @warn "Missing value $x can't be converted to array eltype $T. `missingval` set to `nothing`"
        nothing
    end
end
function _gdalconvertmissing(T::Type{<:Integer}, x::Integer)
    if x >= typemin(T) && x <= typemax(T)
        convert(T, x)
    else
        @warn "Missing value $x can't be converted to array eltype $T. `missingval` set to `nothing`"
        nothing
    end
end
_gdalconvertmissing(T, x) = x

# utility for parsing gdal driver creation options given as Vector{String}
function _parsegdaloptions(options)
    splitted = split.(options, '=')
    splitted = [string.(s) for s in splitted]
    Dict(splitted)
end

function _gdalwrite(filename, A::AbstractRaster, nbands;
    driver=AG.extensiondriver(filename), options=Vector{String}()
)
    A = maybe_typemin_as_missingval(filename, A)
    properties = (width=size(A, X()), height=size(A, Y()), nbands=nbands, dtype=eltype(A))
    if driver isa String
        gdaldriver = AG.getdriver(driver)
    else
        gdaldriver = driver
    end
    options_dict = _parsegdaloptions(options)

    if !("COMPRESS" in keys(options_dict))
        options_dict["COMPRESS"] = "ZSTD"
    end
    
    # drivers supporting the gdal Create() method to directly write to disk
    drivers_supporting_create = ["GTiff", "COG"]  # this should only be a temporary place to put this

    # the goal is to set write block sizes that correspond to eventually blocked reads
    # creation options are driver dependent
    if DA.haschunks(A) == DA.Chunked()
        block_x, block_y = string.(DA.max_chunksize(DA.eachchunk(A)))

        if driver == "GTiff"
            # dont overwrite user specified values
            if !("BLOCKXSIZE" in keys(options_dict))
                options_dict["BLOCKXSIZE"] = block_x
            end
            if !("BLOCKYSIZE" in keys(options_dict))
                options_dict["BLOCKYSIZE"] = block_y
            end
        elseif driver == "COG"
            if !("BLOCKSIZE" in keys(options_dict))
                # cog only supports square blocks
                # if the source already has square blocks, use them
                # otherwise use the driver default
                options_dict["BLOCKSIZE"] = block_x == block_y ? block_x : 512
            end
        end
    end
    # if the input is unchunked we just use the driver defaults

    options_vec = ["$k=$v" for (k,v) in options_dict]

    if driver in drivers_supporting_create
        AG.create(filename; driver=gdaldriver, properties..., options=options_vec) do dataset
            _gdalsetproperties!(dataset, A)
            rds = AG.RasterDataset(dataset)
            open(A; write=true) do O
                rds .= parent(O)
            end
        end
    else
        # Create a memory object and copy it to disk, as ArchGDAL.create
        # does not support direct creation of ASCII etc. rasters
        ArchGDAL.create(""; driver=AG.getdriver("MEM"), properties..., options=options_vec) do dataset
            _gdalsetproperties!(dataset, A)
            rds = AG.RasterDataset(dataset)
            open(A; write=true) do O
                rds .= parent(O)
            end
            AG.copy(dataset; filename=filename, driver=gdaldriver) |> AG.destroy
        end
    end
    return filename
end

function _gdal_bandnames(raster::AG.RasterDataset, nbands = AG.nraster(raster))
    map(1:nbands) do b
        AG.getband(raster.ds, b) do band
            AG.GDAL.gdalgetdescription(band.ptr)
        end
    end
end

function _gdalmetadata(dataset::AG.Dataset, key)
    meta = AG.metadata(dataset)
    regex = Regex("$key=(.*)")
    i = findfirst(f -> occursin(regex, f), meta)
    if i isa Nothing
        return nothing
    else
        return match(regex, meta[i])[1]
    end
end

_gdalsetproperties!(ds, A) = _gdalsetproperties!(ds, dims(A), missingval(A))
function _gdalsetproperties!(dataset, dims, missingval)
    # Convert the dimensions to `Projected` if they are `Converted`
    # This allows saving NetCDF to Tiff
    # Set the index loci to the start of the cell for the lat and lon dimensions.
    # NetCDF or other formats use the center of the interval, so they need conversion.
    x = DD.maybeshiftlocus(GDAL_X_LOCUS, convertlookup(Projected, DD.dims(dims, X)))
    y = DD.maybeshiftlocus(GDAL_Y_LOCUS, convertlookup(Projected, DD.dims(dims, Y)))
    # Convert crs to WKT if it exists
    if !isnothing(crs(x))
        AG.setproj!(dataset, convert(String, convert(WellKnownText, crs(x))))
    end
    # Get the geotransform from the updated lat/lon dims and write
    AG.setgeotransform!(dataset, _dims2geotransform(x, y))

    # Set the nodata value. GDAL can't handle missing. We could choose a default,
    # but we would need to do this for all possible types. `nothing` means
    # there is no missing value.
    # TODO define default nodata values for missing?
    if (!ismissing(missingval) && !isnothing(missingval))
        # We use the axis instead of the values because
        # GDAL has to have values 1:N, not whatever the index holds
        bands = hasdim(dims, Band) ? axes(DD.dims(dims, Band), 1) : 1
        for i in bands
            AG.setnodatavalue!(AG.getband(dataset, i), missingval)
        end
    end

    # Write band labels if they are not Integers.
    if hasdim(dims, Band)
        bandlookup = DD.lookup(dims, Band)
        if !(eltype(bandlookup) <: Integer)
            for i in eachindex(bandlookup)
                AG.getband(dataset, i) do band
                    AG.GDAL.gdalsetdescription(band.ptr, string(bandlookup[i]))
                end
            end
        end
    end

    return dataset
end

# Create a GDAL dataset from any AbstractRaster or dims
function _unsafe_gdal_ds(A::AbstractRaster; missingval=missingval(A), eltype=eltype(A), kw...)
    _unsafe_gdal_ds(dims(A); missingval, eltype, kw...)
end
function _unsafe_gdal_ds(dims::DimTuple; kw...)
    nbands = hasdim(dims, Band) ? length(DD.dims(dims, Band)) : 1
    _unsafe_gdal_ds(dims, nbands; kw...)
end
function _unsafe_gdal_ds(dims::DimTuple, nbands; filename=nothing, suffix=nothing,
    missingval=nothing, metadata=nothing, name=nothing, keys=(name,),
    eltype, driver=_extensiondriver(filename), compress="DEFLATE", chunk=nothing,
)
    gdaldriver = AG.getdriver(driver)
    kw = (
        width=length(DD.dims(dims, X)),
        height=length(DD.dims(dims, Y)),
        nbands=nbands,
        dtype=eltype,
    )
    dataset = if driver == "MEM"
        AG.unsafe_create("tmp"; driver=gdaldriver, kw...)
    elseif driver == "GTiff"
        tileoptions = if isnothing(chunk)
            ["TILED=NO"]
        else
            block_x, block_y = chunk
            ["TILED=YES", "BLOCKXSIZE=$block_x", "BLOCKYSIZE=$block_y"]
        end
        options = ["COMPRESS=$compress", tileoptions...]
        AG.unsafe_create(filename; driver=gdaldriver, options=options, kw...)
    end
    _gdalsetproperties!(dataset, dims, missingval)
    return dataset
end

_extensiondriver(filename::Nothing) = "MEM"
function _extensiondriver(filename::AbstractString)
    # TODO move this check to ArchGDAL
    filename === "/vsimem/tmp" ? "MEM" : AG.extensiondriver(filename)
end

# _maybe_permute_gdal
# Permute dims unless the match the GDAL dimension order
_maybe_permute_to_gdal(A) = _maybe_permute_to_gdal(A, DD.dims(A, (X, Y, Band)))
_maybe_permute_to_gdal(A, dims::Tuple) = A
_maybe_permute_to_gdal(A, dims::Tuple{<:XDim,<:YDim,<:Band}) = permutedims(A, dims)
_maybe_permute_to_gdal(A, dims::Tuple{<:XDim,<:YDim}) = permutedims(A, dims)

_maybe_permute_from_gdal(A, dims::Tuple) = permutedims(A, dims)
_maybe_permute_from_gdal(A, dims::Tuple{<:XDim,<:YDim,<:Band}) = A
_maybe_permute_from_gdal(A, dims::Tuple{<:XDim,<:YDim}) = A

#= Geotranforms ########################################################################

See https://lists.osgeo.org/pipermail/gdal-dev/2011-July/029449.html

"In the particular, but common, case of a “north up” image without any rotation or
shearing, the georeferencing transform takes the following form" :
adfGeoTransform[0] /* top left x */
adfGeoTransform[1] /* w-e pixel resolution */
adfGeoTransform[2] /* 0 */
adfGeoTransform[3] /* top left y */
adfGeoTransform[4] /* 0 */
adfGeoTransform[5] /* n-s pixel resolution (negative value) */
=#

const GDAL_EMPTY_TRANSFORM = [0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
const GDAL_TOPLEFT_X = 1
const GDAL_WE_RES = 2
const GDAL_ROT1 = 3
const GDAL_TOPLEFT_Y = 4
const GDAL_ROT2 = 5
const GDAL_NS_RES = 6

_isalligned(geotransform) = geotransform[GDAL_ROT1] == 0 && geotransform[GDAL_ROT2] == 0

function _geotransform2affine(gt::AbstractVector)
    M = [gt[GDAL_WE_RES] gt[GDAL_ROT1]; gt[GDAL_ROT2] gt[GDAL_NS_RES]]
    v = [gt[GDAL_TOPLEFT_X], gt[GDAL_TOPLEFT_Y]]
    CoordinateTransformations.AffineMap(M, v)
end

function _affine2geotransform(am::CoordinateTransformations.AffineMap)
    M = am.linear
    v = am.translation
    gt = zeros(6)
    gt[GDAL_TOPLEFT_X] = v[1]
    gt[GDAL_WE_RES] = M[1, 1]
    gt[GDAL_ROT1] = M[1, 2]
    gt[GDAL_TOPLEFT_Y] = v[2]
    gt[GDAL_ROT2] = M[2, 1]
    gt[GDAL_NS_RES] = M[2, 2]
    return gt
end

function _dims2geotransform(x::XDim, y::YDim)
    gt = zeros(6)
    gt[GDAL_TOPLEFT_X] = first(x)
    gt[GDAL_WE_RES] = step(x)
    gt[GDAL_ROT1] = zero(eltype(gt))
    gt[GDAL_TOPLEFT_Y] = first(y) - step(y)
    gt[GDAL_ROT2] = zero(eltype(gt))
    gt[GDAL_NS_RES] = step(y)
    return gt
end
function _dims2geotransform(x::XDim{<:AffineProjected}, y::YDim)
    _affine2geotransform(parent(x).affinemap)
end

# precompilation
function _precompile(::Type{GDALfile})
    ccall(:jl_generating_output, Cint, ()) == 1 || return nothing

    for T in (Any, UInt8, UInt16, Int16, UInt32, Int32, Float32, Float64)
        DS = AG.RasterDataset{T,AG.Dataset}
        precompile(crs, (DS,))
        precompile(Rasters.FileArray, (DS, String))
        precompile(dims, (DS,))
        precompile(dims, (DS,WellKnownText{GeoFormatTypes.CRS},Nothing))
        precompile(dims, (DS,WellKnownText{GeoFormatTypes.CRS},EPSG))
        precompile(dims, (DS,WellKnownText{GeoFormatTypes.CRS},ProjString))
        precompile(dims, (DS,WellKnownText{GeoFormatTypes.CRS},WellKnownText{GeoFormatTypes.CRS}))
        precompile(metadata, (DS, key))
        precompile(missingval, (DS, key))
        precompile(Raster, (DS, key))
        precompile(Raster, (DS, String, Nothing))
        precompile(Raster, (DS, String, Symbol))
    end
end

_precompile(GRDfile)
