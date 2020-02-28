using .ArchGDAL

const AG = ArchGDAL

export GDALarray, GDALstack, GDALmetadata, GDALdimMetadata


# Metadata ########################################################################
struct GDALmetadata{K,V} <: ArrayMetadata{K,V}
    val::Dict{K,V}
end

struct GDALdimMetadata{K,V} <: DimMetadata{K,V}
    val::Dict{K,V}
end


# Array ########################################################################

struct GDALarray{T,N,F,D<:Tuple,R<:Tuple,Na<:AbstractString,Me,Mi,W,S
                } <: DiskGeoArray{T,N,D,LazyArray{T,N}}
    filename::F
    dims::D
    refdims::R
    name::Na
    metadata::Me
    missingval::Mi
    window::W
    size::S
end

GDALarray(filename::AbstractString; kwargs...) =
    gdalapply(dataset -> GDALarray(dataset; kwargs...), filename)
GDALarray(dataset::AG.Dataset;
          dims=dims(dataset),
          refdims=(),
          name="",
          metadata=metadata(dataset),
          missingval=missingval(dataset),
          window=()) = begin
    filename = first(AG.filelist(dataset))
    sze = gdalsize(dataset)
    if window != ()
        window = to_indices(dataset, dims2indices(dims, window))
        sze = windowsize(window)
    end
    T = AG.pixeltype(AG.getband(dataset, 1))
    N = length(sze)
    GDALarray{T,N,typeof.((filename,dims,refdims,name,metadata,missingval,window,sze))...
       }(filename, dims, refdims, name, metadata, missingval, window, sze)
end

filename(A::GDALarray) = A.filename
crs(A::GDALarray) = gdalapply(crs, filename(A))
data(A::GDALarray) =
    gdalapply(filename(A)) do dataset
        _window = maybewindow2indices(dataset, dims(A), window(A))
        readwindowed(dataset, _window)
    end

Base.size(A::GDALarray) = A.size

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
    all(hasdim(A, (Lon, Lat))) || error("Array must have Lat and Lon dims to write to GTiff")
    A = permutedims(A, (Lon(), Lat()))
    dataset = AG.unsafe_create(filename;
        width=size(A, 1),
        height=size(A, 2),
        nbands=1,
        dtype=T
    )
    proj = convert(String, crs(dims(A, Lat)))
    AG.setproj!(dataset, proj)
    AG.setgeotransform!(dataset, GDAL_EMPTY_TRANSFORM)
    AG.write!(dataset, data(A), 1)
    AG.destroy(dataset)
    return filename
end
Base.write(filename::AbstractString, ::Type{GDALarray}, A::GeoArray{T,3}) where T = begin
    DimensionalData.hasdim(A, Band()) || error("Must have a `Band` dimension to write a 3-dimensional array")
    nbands = size(A, Band())
    A = permutedims(A, (Lon(), Lat(), Band()))
    dataset = AG.unsafe_create(filename;
        width=size(A, 1),
        height=size(A, 2),
        nbands=nbands,
        dtype=T,
    )
    proj = convert(String, crs(dims(A, Lat)))
    AG.setgeotransform!(dataset, GDAL_EMPTY_TRANSFORM)
    AG.setproj!(dataset, proj)
    AG.write!(dataset, data(A), Cint[1])
    AG.destroy(dataset)
    return filename
end


# Stack ########################################################################

struct GDALstack{T,R,W,M} <: DiskGeoStack{T}
    filename::T
    refdims::R
    window::W
    metadata::M
end

GDALstack(filenames::NamedTuple;
          refdims=(), window=(),
          metadata=gdalapply(metadata, first(values(filenames)))) =
    GDALstack(filenames, refdims, window, metadata)

safeapply(f, ::GDALstack, path::AbstractString) = gdalapply(f, path)

@inline Base.getindex(s::GDALstack, key::Key) =
    gdalapply(filename(s, key)) do dataset
        GDALarray(dataset; refdims=refdims(s), name=string(key), window=window(s))
    end
@inline Base.getindex(s::GDALstack, key::Key, I::Union{Colon,Integer,AbstractArray}...) =
    s[key][I...]


Base.copy!(dst::AbstractGeoArray, src::GDALstack, key::Key) =
    copy!(data(dst), src, key)
Base.copy!(dst::AbstractArray, src::GDALstack, key::Key) =
    gdalapply(filename(src, key)) do dataset
        key = string(key)
        _window = maybewindow2indices(dataset, dims(dataset), window(src))
        copy!(dst, readwindowed(dataset, _window))
    end


# DimensionalData methods for ArchGDAL types ###############################

dims(dataset::AG.Dataset) = begin
    gt = try
        AG.getgeotransform(dataset)
    catch
        GDAL_EMPTY_TRANSFORM
    end

    ysize, xsize = AG.height(dataset), AG.width(dataset)

    nbands = AG.nraster(dataset)
    band = Band(1:nbands, grid=CategoricalGrid())
    sourcecrs = crs(dataset)
    targetcrs = EPSG(4326)

    lonlat_metadata=Dict(:crs => sourcecrs)

    # Output a BoundedGrid dims when the transformation is lat/lon alligned,
    # otherwise use TransformedGrid with an affine map.
    if isalligned(gt)
        lonspan = lonres(gt)

        lonmin = gt[GDAL_TOPLEFT_X]
        lonmax = lonmin + lonspan * xsize
        loncoords = reproject(tuple.(LinRange(lonmin, lonmax - lonspan, xsize), 0.0), sourcecrs, targetcrs)
        if !isapprox(loncoords[1][1], loncoords[end][1]; rtol=1e-5)
            error("Longitude dimension is not grid-alligned $(loncoords[1][1]) $(loncoords[end][1])")
        end
        lonrange = last.(loncoords)
        lonbounds = lonrange[1], reproject([(lonmax, 0.0)], sourcecrs, targetcrs)[1][1]
        # lonbounds = lonmin, lonmax
        lon = Lon(lonrange; grid=BoundedGrid(bounds=lonbounds), metadata=lonlat_metadata)

        latspan = latres(gt)
        latmax = gt[GDAL_TOPLEFT_Y]
        latmin = latmax + latspan * (ysize - 1)
        latcoords = reproject(tuple.(0.0, LinRange(latmin, latmax, ysize)), sourcecrs, targetcrs)
        if !isapprox(latcoords[1][2], latcoords[end][2]; rtol=1e-5)
            error("Latitude dimension is not grid-alligned $(latcoords[1][2]) $(latcoords[end][2])")
        end
        latrange = first.(latcoords)
        latbounds = latrange[1], reproject([(0.0, latmax)], sourcecrs, targetcrs)[1][2]
        # latbounds = latmin, latmax
        latgrid = BoundedGrid(order=Ordered(Forward(), Reverse(), Reverse()), bounds=latbounds)
        lat = Lat(latrange; grid=latgrid, metadata=lonlat_metadata)

        formatdims((1:xsize, 1:ysize, 1:nbands), (lon, lat, band))
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
    T = AG.pixeltype(band)
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

crs(dataset::AG.Dataset, args...) =
    WellKnownText(GeoFormatTypes.CRS(), string(AG.getproj(dataset)))

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
    AG.read(filepath) do dataset
        f(dataset)
    end

gdalread(s::GDALstack, key, I...) =
    gdalapply(filename(s, key)) do dataset
        readwindowed(dataset, window(s), I...)
    end
gdalread(A::GDALarray, I...) =
    gdalapply(filename(A)) do dataset
        readwindowed(dataset, window(A), I...)
    end

gdalsize(dataset) = begin
    band = AG.getband(dataset, 1)
    AG.width(band), AG.height(band), AG.nraster(dataset)
end

# See https://lists.osgeo.org/pipermail/gdal-dev/2011-July/029449.html
# for an explanation of the geotransform format
geotransform_to_affine(gt) = begin
    AffineMap([gt[2] gt[3]; gt[5] gt[6]], [gt[1], gt[4]])
end

reproject(coords, crs, taaaarget) =  begin
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
