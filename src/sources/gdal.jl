using ArchGDAL

const AG = ArchGDAL

export GDALarray, GDALstack, GDALmetadata, GDALdimMetadata


# Metadata ########################################################################

struct GDALmetadata{M} <: AbstractArrayMetadata
    val::M
end

struct GDALdimMetadata{M} <: AbstractDimMetadata
    val::M
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
          selector_crs=nothing,
          dims=dims(dataset; selector_crs=selector_crs),
          refdims=(),
          metadata=metadata(dataset),
          missingval=missingval(dataset),
          name="Unnamed",
          window=(),
         ) = begin
    filename = first(AG.filelist(dataset))
    if window == ()
        sze = gdalsize(dataset)
    else
        sze = windowsize(window)
        dims, refdims = slicedims(dims, refdims, window)
    end
    T = AG.getdatatype(AG.getband(dataset, 1))
    N = length(sze)
    try
        missingval = convert(T, missingval)
    catch
        @warn "No data value from GDAL $(missingval) is not convertible to data type $T. `missingval` is probably incorrect."
    end
    GDALarray{T,N,typeof.((filename,dims,refdims,metadata,missingval,name,window,sze))...
       }(filename, dims, refdims, metadata, missingval, name, window, sze)
end

Base.size(A::GDALarray) = A.size

Base.parent(A::GDALarray) =
    gdalapply(dataset -> gdalread(dataset, windoworempty(A)...), A.filename)

Base.getindex(A::GDALarray, I::Vararg{<:Union{<:Integer,<:AbstractArray}}) = begin
    I = applywindow(A, I)
    rebuildsliced(A, gdalapply(dataset -> gdalread(dataset, I...), A.filename), I)
end


# Stack ########################################################################

struct GDALstack{T,D,R,W,M} <: AbstractGeoStack{T}
    data::T
    dims::D
    refdims::R
    window::W
    metadata::M
end

GDALstack(data::NamedTuple;
          dims=gdalapply(dims, first(values(data))),
          refdims=(), window=(),
          metadata=gdalapply(metadata, first(values(data)))) =
    GDALstack(data, dims, refdims, window, metadata)


@inline rebuild(s::GDALstack; data=parent(s), dims=dims(s), refdims=refdims(s),
        window=window(s), metadata=metadata(s)) =
    GDALstack(data, dims, refdims, window, metadata)

safeapply(f, ::GDALstack, path::AbstractString) = gdalapply(f, path)

data(::GDALstack, dataset, key::Key, I...) = GDALarray(dataset; window=I)
data(::GDALstack, dataset, key::Key) = GDALarray(dataset)


Base.write(filename::AbstractString, ::Type{<:AbstractGeoArray}, a::AbstractGeoArray) =
    Base.write(GDALarray, filename, GeoArray(a))
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
        proj = convert(String, crs(dims(A, Lon)))
        AG.setproj!(dataset, proj)
        AG.setgeotransform!(dataset, GDAL_EMPTY_TRANSFORM)
        AG.write!(dataset, source(A), 1)
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
        crs = convert(String, crs(A, Lon))
        AG.setgeotransform!(dataset, GDAL_EMPTY_TRANSFORM)
        AG.setproj!(dataset, proj)
        AG.write!(dataset, source(A), Cint[1])
        AG.destroy(dataset)
    end
end
Base.write(filename::AbstractString, ::Type{GDALarray}, s::AbstractGeoStack) =
    for key in keys(s)
        fn = joinpath(dirname(filename), string(key, "_", basename(filename)))
        write(fn, GDALarray, s[key])
    end

Base.copy!(dst::AbstractArray, src::GDALstack, key::Key) =
    copy!(dst, gdalapply(AG.read, source(src, key)))
Base.copy!(dst::AbstractGeoArray, src::GDALstack, key::Key) =
    copy!(parent(dst), gdalapply(AG.read, source(src, key)))



# DimensionalData methods for ArchGDAL types ###############################

dims(dataset::AG.Dataset; selector_crs=nothing) = begin
    gt = AG.getgeotransform(dataset)
    # gt = GDAL_EMPTY_TRANSFORM
    ysize, xsize = AG.height(dataset), AG.width(dataset)
    crs = WellKnownText(string(AG.getproj(dataset)))

    nbands = AG.nraster(dataset)
    band = Band(1:nbands, grid=CategoricalGrid())

    # Output a AllignedGrid dims when the transformation is lat/lon alligned,
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
        lon = Lon(lonrange; 
                  grid=RegularGrid(span=abs(lonspan), crs=crs, selector_crs=selector_crs))

        latspan = latres(gt)
        latmax = gt[GDAL_TOPLEFT_Y]
        latmin = latmax + latspan * (ysize - 1)
        latcoords = reproject(tuple.(0.0, LinRange(latmin, latmax, ysize)), crs)
        if latcoords[1][1] != latcoords[end][1]
            error("Latitude dimension is not grid-alligned $(latcoords[1][1]) $(latcoords[end][1])")
        end
        latrange = last.(latcoords)
        lat = Lat(latrange;
                  grid=RegularGrid(order=Ordered(Forward(), Reverse(), Forward()),
                                   span=abs(latspan), crs=crs, selector_crs=selector_crs))

        formatdims((xsize, ysize, nbands), (lon, lat, band))
    else
        error("Rotated grids not handled currently")
        # affinemap = geotransform_to_affine(geotransform)
        # x = X(affinemap; grid=TransformedGrid(dims=Lon()))
        # y = Y(affinemap; grid=TransformedGrid(dims=Lat()))

        # formatdims((xsize, ysize, nbands), (x, y, band))
    end
end

missingval(dataset::AG.Dataset, args...) =
    AG.getnodatavalue(AG.getband(dataset, 1))
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


gdalapply(f, path::AbstractString) =
    AG.registerdrivers() do
        AG.read(path) do dataset::AG.Dataset
            f(dataset)
        end
    end

gdalread(dataset, I...) = AG.read(dataset, reverse(I)...)
gdalread(dataset) = AG.read(dataset)

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
