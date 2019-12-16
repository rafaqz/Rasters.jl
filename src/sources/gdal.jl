using ArchGDAL

const AG = ArchGDAL

export GDALarray, GDALstack, GDALmetadata, GDALdimMetadata


# Metadata ########################################################################

struct GDALmetadata{K,V} <: AbstractArrayMetadata{K,V}
    val::Dict{K,V}
end

struct GDALdimMetadata{K,V} <: AbstractDimMetadata{K,V}
    val::Dict{K,V}
end


# Array ########################################################################

struct GDALarray{T,N,A,D<:Tuple,R<:Tuple,Me,Mi,Na,W,S} <: DiskGeoArray{T,N,D}
    filename::A
    dims::D
    refdims::R
    metadata::Me
    missingval::Mi
    name::Na
    window::W
    size::S
end

GDALarray(filename::AbstractString; kwargs...) =
    gdalapply(dataset -> GDALarray(dataset; kwargs...), filename)
GDALarray(dataset::AG.Dataset;
          dims=dims(dataset),
          refdims=(),
          metadata=metadata(dataset),
          missingval=missingval(dataset),
          name="Unnamed",
          window=()) = begin
    filename = first(AG.filelist(dataset))
    if window == ()
        sze = gdalsize(dataset)
    else
        window = dims2indices(dims, window)
        sze = windowsize(window)
    end
    T = AG.getdatatype(AG.getband(dataset, 1))
    N = length(sze)
    GDALarray{T,N,typeof.((filename,dims,refdims,metadata,missingval,name,window,sze))...
       }(filename, dims, refdims, metadata, missingval, name, window, sze)
end

Base.size(A::GDALarray) = A.size


Base.parent(A::GDALarray) =
    gdalapply(filename(A)) do dataset
        _window = maybewindow2indices(dataset, dims(A), window(A))
        readwindowed(dataset, _window)
    end
Base.getindex(A::GDALarray, I::Vararg{<:Union{<:Integer,<:AbstractArray}}) =
    gdalapply(filename(A)) do dataset
        _window = maybewindow2indices(dataset, dims(A), window(A))
        # Slice for both window and indices
        _dims, _refdims = slicedims(slicedims(dims(A), refdims(A), _window)..., I)
        data = readwindowed(dataset, _window, I...)
        rebuild(A, data, _dims, _refdims)
    end
Base.getindex(A::GDALarray, i1::Integer, I::Vararg{<:Integer}) =
    gdalapply(filename(A)) do dataset
        _window = maybewindow2indices(dataset, dims(A), window(A))
        readwindowed(dataset, _window, i1, I...)
    end

Base.write(filename::AbstractString, ::Type{GDALarray}, A::GeoArray{T,2}) where T = begin
    AG.registerdrivers() do
        driver = AG.getdriver("GTiff")
        A = permutedims(A, (Lon(), Lat()))
        dataset = AG.unsafe_create(filename, driver;
            width = size(A, 1),
            height = size(A, 2),
            nbands = 1,
            dtype = T
        )
        proj = "" #convert(String, crs(A))
        AG.setproj!(dataset, proj)
        AG.setgeotransform!(dataset, GDAL_EMPTY_TRANSFORM)
        AG.write!(dataset, parent(A), 1)
        AG.destroy(dataset)
    end
end
Base.write(filename::AbstractString, ::Type{GDALarray}, A::GeoArray{T,3}) where T = begin
    DimensionalData.hasdim(A, Band()) || error("Must have a `Band` dimension to write a 3-dimensional array")
    nbands = size(A, Band())
    AG.registerdrivers() do
        driver = AG.getdriver("GTiff") # Returns NULL Driver?
        A = permutedims(A, (Lon(), Lat(), Band()))
        dataset = AG.unsafe_create(filename, driver;
            width = size(A, 1),
            height = size(A, 2),
            nbands = nbands,
            dtype = T
        )
        proj = convert(String, projection(A))
        AG.setgeotransform!(dataset, GDAL_EMPTY_TRANSFORM)
        AG.setproj!(dataset, proj)
        AG.write!(dataset, parent(A), Cint[1])
        AG.destroy(dataset)
    end
end


# Stack ########################################################################

struct GDALstack{T,D,R,W,M} <: DiskGeoStack{T}
    filename::T
    dims::D
    refdims::R
    window::W
    metadata::M
end

GDALstack(filenames::NamedTuple;
          dims=gdalapply(dims, first(values(filenames))),
          refdims=(), window=(),
          metadata=gdalapply(metadata, first(values(filenames)))) =
    GDALstack(filenames, dims, refdims, window, metadata)


@inline rebuild(s::GDALstack; data=filename(s), dims=dims(s), refdims=refdims(s),
        window=window(s), metadata=metadata(s)) =
    GDALstack(data, dims, refdims, window, metadata)

safeapply(f, ::GDALstack, path::AbstractString) = gdalapply(f, path)

@inline Base.getindex(s::GDALstack, key::Key, i1::Integer, I::Integer...) =
    gdalapply(filename(s, key)) do dataset
        _window = maybewindow2indices(dataset, dims(s), window(s))
        readwindowed(dataset, _window, I...)
    end
@inline Base.getindex(s::GDALstack, key::Key, I::Union{Colon,Integer,AbstractArray}...) =
    gdalapply(filename(s, key)) do dataset
        _dims = dims(s)
        _window = maybewindow2indices(dataset, _dims, window(s))
        _dims, _refdims = slicedims(slicedims(_dims, refdims(s), _window)..., I)
        A = readwindowed(dataset, _window, I...)
        GeoArray(A, _dims, _refdims, metadata(dataset), missingval(dataset), string(key))
    end


Base.copy!(dst::AbstractGeoArray, src::GDALstack, key::Key) =
    copy!(parent(dst), src, key)
Base.copy!(dst::AbstractArray, src::GDALstack, key::Key) =
    gdalapply(filename(src, key)) do dataset
        key = string(key)
        _window = maybewindow2indices(dataset, dims(dataset), window(src))
        copy!(dst, readwindowed(dataset, _window))
    end


# DimensionalData methods for ArchGDAL types ###############################

dims(dataset::AG.Dataset) = begin
    gt = AG.getgeotransform(dataset)
    # gt = GDAL_EMPTY_TRANSFORM
    ysize, xsize = AG.height(dataset), AG.width(dataset)
    crs = WellKnownText(string(AG.getproj(dataset)))

    nbands = AG.nraster(dataset)
    band = Band(1:nbands, grid=CategoricalGrid())

    # Output a BoundedGrid dims when the transformation is lat/lon alligned,
    # otherwise use TransformedGrid with an affine map.
    if isalligned(gt)
        lonspan = lonres(gt)

        lonmin = gt[GDAL_TOPLEFT_X]
        lonmax = lonmin + lonspan * xsize
        loncoords = reproject(tuple.(LinRange(lonmin, lonmax - lonspan, xsize), 0.0), crs)
        if loncoords[1][2] != loncoords[end][2]
            error("Longitude dimension is not grid-alligned $(loncoords[1][2]) $(loncoords[end][2])")
        end
        lonrange = first.(loncoords)
        lonbounds = lonrange[1], reproject([(lonmax, 0.0)], crs)[1][1]
        lon = Lon(lonrange; grid=BoundedGrid(bounds=lonbounds))

        latspan = latres(gt)
        latmax = gt[GDAL_TOPLEFT_Y]
        latmin = latmax + latspan * (ysize - 1)
        latcoords = reproject(tuple.(0.0, LinRange(latmin, latmax, ysize)), crs)
        if latcoords[1][1] != latcoords[end][1]
            error("Latitude dimension is not grid-alligned $(latcoords[1][1]) $(latcoords[end][1])")
        end
        latrange = last.(latcoords)
        latbounds = latrange[1], reproject([(0.0, latmax)], crs)[1][2]
        latgrid = BoundedGrid(order=Ordered(Forward(), Reverse(), Forward()), bounds=latbounds)
        lat = Lat(latrange; grid=latgrid)

        formatdims((xsize, ysize, nbands), (lon, lat, band))
    else
        error("Rotated grids not handled currently")
        # affinemap = geotransform_to_affine(geotransform)
        # x = X(affinemap; grid=TransformedGrid(dims=Lon()))
        # y = Y(affinemap; grid=TransformedGrid(dims=Lat()))

        # formatdims((xsize, ysize, nbands), (x, y, band))
    end
end

missingval(dataset::AG.Dataset, args...) = begin
    band = AG.getband(dataset, 1)
    missingval = AG.getnodatavalue(band)
    T = AG.getdatatype(band)
    try
        missingval = convert(T, missingval)
    catch
        @warn "No data value from GDAL $(missingval) is not convertible to data type $T. `missingval` is probably incorrect."
    end
    missingval
end
metadata(dataset::AG.Dataset, args...) = begin
    band = AG.getband(dataset, 1)
    # color = AG.getname(AG.getcolorinterp(band))
    scale = AG.getscale(band)
    offset = AG.getoffset(band)
    # norvw = AG.noverview(band)
    units = AG.getunittype(band)
    path = first(AG.filelist(dataset))
    GDALmetadata(Dict("filepath"=>path, "scale"=>scale, "offset"=>offset, "units"=>units))
end


# Utils ########################################################################

#=
In the particular, but common, case of a “north up” image without any rotation or shearing, the georeferencing transform takes the following form :
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

isalligned(geotransform) = geotransform[GDAL_ROT1] == 0 && geotransform[GDAL_ROT2] == 0
latres(geotransform) = geotransform[GDAL_NS_RES]
lonres(geotransform) = geotransform[GDAL_WE_RES]


gdalapply(f, filepath::AbstractString) =
    AG.registerdrivers() do
        AG.read(filepath) do dataset
            f(dataset)
        end
    end

gdalread(s::GDALstack, key, I...) =
    gdalapply(filename(s, key)) do dataset
        readwindowed(dataset, window(s), I...)
    end
gdalread(A::GDALarray, I...) =
    gdalapply(filename(A)) do dataset
        readwindowed(dataset, window(A), I...)
    end

readwindowed(A::AG.Dataset, window::Tuple{}) = AG.read(A)

gdalsize(dataset) = begin
    band = AG.getband(dataset, 1)
    AG.width(band), AG.height(band), AG.nraster(dataset)
end

# See https://lists.osgeo.org/pipermail/gdal-dev/2011-July/029449.html
# for an explanation of the geotransform format
geotransform_to_affine(gt) = begin
    AffineMap([gt[2] gt[3]; gt[5] gt[6]], [gt[1], gt[4]])
end

reproject(coords, crs) =  begin
    AG.importWKT(GeoFormatTypes.val(crs)) do source
        AG.importEPSG(4326) do target
            AG.createcoordtrans(source, target) do transform
                transformcoord.(coords, Ref(transform))
            end
        end
    end
end

transformcoord(coord, transform) =
    AG.createpoint(coord...) do point
        AG.coordinates(AG.transform!(point, transform))
    end
