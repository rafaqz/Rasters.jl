const AG = ArchGDAL

const GDAL_X_INDEX = ForwardIndex()
const GDAL_Y_INDEX = ReverseIndex()
const GDAL_BAND_INDEX = ForwardIndex()
const GDAL_X_ARRAY = ForwardArray()
const GDAL_Y_ARRAY = ReverseArray()
const GDAL_BAND_ARRAY = ForwardArray()
const GDAL_RELATION = ForwardRelation()
const GDAL_X_LOCUS = Start()
const GDAL_Y_LOCUS = Start()

# Array ########################################################################

FileArray{_GRD}(filename) = _read(ds -> FileArray(ds, filename), _GDAL, filename)
function FileArray(raster::AG.RasterDataset, filename, key=nothing)
    FileArray{_GDAL,eltype(raster),ndims(raster)}(filename, size(raster))
end

Base.open(f::Function, A::FileArray{_GDAL}, key...) = _read(f, _GDAL, filename(A))

# AbstractGeoArray methods

"""
    Base.write(filename::AbstractString, ::Type{GDALarray}, A::AbstractGeoArray; kw...)

Write a [`GDALarray`](@ref) to file, `.tif` by default, but other GDAL drivers also work.

# Keywords

- `driver::String`: a GDAL driver name. Guessed from the filename extension by default.
- `compress::String`: GeoTIFF compression flag. "DEFLATE" by default.
- `tiled::Bool`: GeoTiff tiling. Defaults to `true`.

Returns `filename`.
"""
function Base.write(
    filename::AbstractString, ::Type{_GDAL}, A::AbstractGeoArray{T,2}; kw...
) where T
    all(hasdim(A, (XDim, Y))) || error("Array must have Y and X dims")

    correctedA = permutedims(A, (X(), Y())) |>
        a -> reorder(a, (X(GDAL_X_INDEX), Y(GDAL_Y_INDEX))) |>
        a -> reorder(a, GDAL_RELATION)
    checkarrayorder(correctedA, (GDAL_X_ARRAY, GDAL_Y_ARRAY))
    checkindexorder(correctedA, (GDAL_X_INDEX, GDAL_Y_INDEX))

    nbands = 1
    indices = 1
    _gdalwrite(filename, correctedA, nbands, indices; kw...)
end
function Base.write(
    filename::AbstractString, ::Type{_GDAL}, A::AbstractGeoArray{T,3}, kw...
) where T
    all(hasdim(A, (X, Y))) || error("Array must have Y and X dims")
    hasdim(A, Band()) || error("Must have a `Band` dimension to write a 3-dimensional array")

    correctedA = permutedims(A, (X(), Y(), Band())) |>
        a -> reorder(a, (X(GDAL_X_INDEX), Y(GDAL_Y_INDEX), Band(GDAL_BAND_INDEX))) |>
        a -> reorder(a, GDAL_RELATION)
    checkarrayorder(correctedA, (GDAL_X_ARRAY, GDAL_Y_ARRAY, GDAL_BAND_ARRAY))
    checkindexorder(correctedA, (GDAL_X_INDEX, GDAL_Y_INDEX, GDAL_BAND_INDEX))

    nbands = size(correctedA, Band())
    indices = Cint[1:nbands...]
    _gdalwrite(filename, correctedA, nbands, indices; kw...)
end


# DimensionalData methods for ArchGDAL types ###############################

function DD.dims(raster::AG.RasterDataset, crs=nothing, mappedcrs=nothing)
    gt = try
        AG.getgeotransform(raster) catch GDAL_EMPTY_TRANSFORM end
    xsize, ysize = size(raster)

    nbands = AG.nraster(raster)
    band = Band(1:nbands, mode=Categorical(Ordered()))
    crs = crs isa Nothing ? GeoData.crs(raster) : crs
    xy_metadata = Metadata{_GDAL}()

    # Output Sampled index dims when the transformation is lat/lon alligned,
    # otherwise use Transformed index, with an affine map.
    if _isalligned(gt)
        xstep = gt[GDAL_WE_RES]
        xmin = gt[GDAL_TOPLEFT_X]
        xmax = gt[GDAL_TOPLEFT_X] + xstep * (xsize - 1)
        xindex = LinRange(xmin, xmax, xsize)

        ystep = gt[GDAL_NS_RES] # A negative number
        ymax = gt[GDAL_TOPLEFT_Y] + ystep
        ymin = gt[GDAL_TOPLEFT_Y] + ystep * ysize
        yindex = LinRange(ymax, ymin, ysize)

        # Spatial data defaults to area/inteval
        xsampling, ysampling = if _gdalmetadata(raster.ds, "AREA_OR_POINT") == "Point"
            Points(), Points()
        else
            # GeoTiff uses the "pixelCorner" convention
            Intervals(GDAL_X_LOCUS), Intervals(GDAL_Y_LOCUS)
        end

        xmode = Projected(
            order=Ordered(GDAL_X_INDEX, GDAL_X_ARRAY, GDAL_RELATION),
            span=Regular(step(xindex)),
            sampling=xsampling,
            crs=crs,
            mappedcrs=mappedcrs,
        )
        ymode = Projected(
            order=Ordered(GDAL_Y_INDEX, GDAL_Y_ARRAY, GDAL_RELATION),
            sampling=ysampling,
            # Use the range step as is will be different to ystep due to float error
            span=Regular(step(yindex)),
            crs=crs,
            mappedcrs=mappedcrs,
        )
        x = X(xindex; mode=xmode, metadata=xy_metadata)
        y = Y(yindex; mode=ymode, metadata=xy_metadata)

        DimensionalData._formatdims(map(Base.OneTo, (xsize, ysize, nbands)), (x, y, band))
    else
        error("Rotated/transformed dimensions are not handled yet. Open a github issue for GeoData.jl if you need this.")
        # affinemap = geotransform2affine(geotransform)
        # x = X(affinemap; mode=TransformedIndex(dims=X()))
        # y = Y(affinemap; mode=TransformedIndex(dims=Y()))

        # formatdims((xsize, ysize, nbands), (x, y, band))
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
    upair = units == "" ? () : (:units=>units,)
    Metadata{_GDAL}(Dict(:filepath=>path, :scale=>scale, :offset=>offset, upair...))
end

function missingval(raster::AG.RasterDataset, args...)
    # We can only handle data where all bands have the same missingval
    band = AG.getband(raster.ds, 1)
    nodata = AG.getnodatavalue(band)
    return if nodata isa Nothing
        nothing
    else
        # Convert in case getnodatavalue is the wrong type.
        # This conversion should always be safe.
        if eltype(band) <: AbstractFloat && nodata isa Real
            convert(eltype(band), nodata)
        else
            nodata
        end
    end
end

crs(raster::AG.RasterDataset, args...) =
    WellKnownText(GeoFormatTypes.CRS(), string(AG.getproj(raster.ds)))


# Utils ########################################################################

function _gdalmetadata(dataset::AG.Dataset, key)
    meta = AG.metadata(dataset)
    regex = Regex("$key=(.*)")
    i = findfirst(f -> occursin(regex, f), meta)
    if i isa Nothing
        nothing
    else
        match(regex, meta[i])[1]
    end
end

function _read(f, ::Type{_GDAL}, filename::AbstractString, key...)
    AG.readraster(filename) do raster
        f(raster)
    end
end

function _gdalwrite(filename, A, nbands, indices; 
    driver=AG.extensiondriver(filename), compress="DEFLATE", tiled=true
)
    tiledstring = tiled isa Bool ? (tiled ? "YES" : "NO") : tiled
    kw = (width=size(A, X()), height=size(A, Y()), nbands=nbands, dtype=eltype(A))
    gdaldriver = AG.getdriver(driver)
    if driver == "GTiff" 
        options = ["COMPRESS=$compress", "TILED=$tiledstring"]
        AG.create(filename; driver=gdaldriver, options=options, kw...) do dataset
            _gdalsetproperties!(dataset, A)
            AG.write!(dataset, readwindowed(A), indices)
        end
    else
        # Create a  memory object and copy it to disk, as ArchGDAL.create
        # does not support direct creation of ASCII etc. rasters
        ArchGDAL.create(""; driver=AG.getdriver("MEM"), kw...) do dataset
            _gdalsetproperties!(dataset, A)
            AG.write!(dataset, readwindowed(A), indices)
            AG.copy(dataset; filename=filename, driver=gdaldriver) |> AG.destroy
        end
    end
    return filename
end

function _gdalsetproperties!(dataset, A)
    # Convert the dimensions to `Projected` if they are `Converted`
    # This allows saving NetCDF to Tiff
    # Set the index loci to the start of the cell for the lat and lon dimensions.
    # NetCDF or other formats use the center of the interval, so they need conversion.
    x = convertmode(Projected, DD.maybeshiftlocus(GDAL_X_LOCUS, dims(A, X)))
    y = convertmode(Projected, DD.maybeshiftlocus(GDAL_Y_LOCUS, dims(A, Y)))
    # Convert crs to WKT if it exists
    if !(crs(x) isa Nothing)
        AG.setproj!(dataset, convert(String, convert(WellKnownText, crs(x))))
    end
    # Get the geotransform from the updated lat/lon dims and write
    AG.setgeotransform!(dataset, _dims2geotransform(x, y))

    # Set the nodata value. GDAL can't handle missing. We could choose a default, 
    # but we would need to do this for all possible types. `nothing` means
    # there is not missing value.
    # TODO define default nodata values for missing
    if (missingval(A) !== missing) && (missingval(A) !== nothing)
        bands = hasdim(A, Band) ? index(A, Band) : 1
        for i in bands
            AG.setnodatavalue!(AG.getband(dataset, i), missingval(A))
        end
    end

    dataset
end

# Create a GeoArray from a memory-backed dataset
function GeoArray(dataset::AG.Dataset;
    crs=nothing, mappedcrs=nothing,
    dims=dims(AG.RasterDataset(dataset), crs, mappedcrs),
    refdims=(), name=Symbol(""),
    metadata=metadata(AG.RasterDataset(dataset)),
    missingval=missingval(AG.RasterDataset(dataset))
)
    GeoArray(AG.read(dataset), dims, refdims, name, metadata, missingval)
end

# Create a memory-backed GDAL dataset from any AbstractGeoArray
function unsafe_ArchGDALdataset(A::AbstractGeoArray)
    width = size(A, X)
    height = size(A, Y)
    nbands = size(A, Band)

    dataset = AG.unsafe_create("tmp";
        driver=AG.getdriver("MEM"),
        width=width,
        height=height,
        nbands=nbands,
        dtype=eltype(A)
    )
    # write bands to dataset
    AG.write!(dataset, data(permutedims(A, (X, Y, Band))), Cint[1:nbands...])
    _gdalsetproperties!(dataset, A)
    dataset
end

function AG.Dataset(f::Function, A::AbstractGeoArray)
    dataset = unsafe_ArchGDALdataset(A)
    try
        f(dataset)
    finally
        AG.destroy(dataset)
    end
end


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

_geotransform2affine(gt) =
    AffineMap([gt[GDAL_WE_RES] gt[GDAL_ROT1]; gt[GDAL_ROT2] gt[GDAL_NS_RES]],
              [gt[GDAL_TOPLEFT_X], gt[GDAL_TOPLEFT_Y]])

function _dims2geotransform(x::X, y::Y)
    gt = zeros(6)
    gt[GDAL_TOPLEFT_X] = first(x)
    gt[GDAL_WE_RES] = step(x)
    gt[GDAL_ROT1] = 0.0
    gt[GDAL_TOPLEFT_Y] = first(y) - step(y)
    gt[GDAL_ROT2] = 0.0
    gt[GDAL_NS_RES] = step(y)
    return gt
end
