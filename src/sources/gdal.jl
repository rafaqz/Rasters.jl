using .ArchGDAL

const AG = ArchGDAL

export GDALarray, GDALstack

# Array ########################################################################

@GeoArrayMixin struct GDALarray{W,S,A} <: AbstractGeoArray{T,N,D}
    window::W
    size::S
end

GDALarray(path::AbstractString; kwargs...) =
    gdalapply(dataset -> GDALarray(dataset; kwargs...), path)
GDALarray(dataset::AG.Dataset;
          dims=dims(dataset),
          refdims=(),
          metadata=metadata(dataset),
          missingval=missingval(dataset),
          name=Symbol(""),
          window=()) = begin
    path = first(AG.filelist(dataset))
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
    GDALarray{T,N,typeof.((dims,refdims,metadata,missingval,name,window,sze,path))...
       }(path, dims, refdims, metadata, missingval, name, window, sze)
end

Base.size(A::GDALarray) = A.size
Base.parent(A::GDALarray) =
    gdalapply(dataset -> gdalread(dataset, windoworempty(A)...), A.data)
Base.getindex(A::GDALarray, I::Vararg{<:Union{<:Integer,<:AbstractArray}}) = begin
    I = applywindow(A, I)
    rebuildsliced(A, gdalapply(dataset -> gdalread(dataset, I...), A.data), I)
end

# Stack ########################################################################

@GeoStackMixin struct GDALstack{} <: AbstractGeoStack{T} end

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

Base.copy!(dst::AbstractArray, src::GDALstack, key::Key) =
    copy!(dst, gdalapply(AG.read, source(src, key)))
Base.copy!(dst::AbstractGeoArray, src::GDALstack, key::Key) =
copy!(parent(dst), gdalapply(AG.read, source(src, key)))


# DimensionalData methods for ArchGDAL types ###############################

dims(dataset::AG.Dataset) = begin
    gt = AG.getgeotransform(dataset)
    ysize, xsize = AG.height(dataset), AG.width(dataset)
    crs = WellKnownText(string(AG.getproj(dataset)))

    nbands = AG.nraster(dataset)
    band = Band(1:nbands, grid=CategoricalGrid())

    # Output a AllignedGrid dims when the transformation is lat/lon alligned,
    # otherwise use TransformedGrid with an affine map.
    if isalligned(gt)
        lonspan = lonres(gt)

        lonmin = gt[GDAL_TOPLEFT_X]
        lonmax = lonmin + lonspan * (xsize - 1)
        loncoords = reproject(tuple.(LinRange(lonmin, lonmax, xsize), 0.0), crs)
        if loncoords[1][2] != loncoords[end][2]
            error("Longitude dimension is not grid-alligned $(loncoords[1][2]) $(loncoords[end][2])")
        end
        lonrange = first.(loncoords)
        longrid = AllignedGrid(span=abs(lonspan))
        lon = Lon(lonrange; grid=longrid)

        latspan = latres(gt)
        latmax = gt[GDAL_TOPLEFT_Y]
        latmin = latmax + latspan * (ysize - 1)
        latcoords = reproject(tuple.(0.0, LinRange(latmin, latmax, ysize)), crs)
        if latcoords[1][1] != latcoords[end][1] 
            error("Latitude dimension is not grid-alligned $(latcoords[1][1]) $(latcoords[end][1])")
        end
        latrange = last.(latcoords) 
        latgrid = AllignedGrid(order=Ordered(Forward(), Reverse()), span=abs(latspan))
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

missingval(dataset::AG.Dataset, args...) =
    AG.getnodatavalue(AG.getband(dataset, 1))
metadata(dataset::AG.Dataset, args...) = begin
    band = AG.getband(dataset, 1)
    color = AG.getname(AG.getcolorinterp(band))
    scale = AG.getscale(band)
    offset = AG.getoffset(band)
    norvw = AG.noverview(band)
    units = AG.getunittype(band)
    crs = WellKnownText(string(AG.getproj(dataset)))
    path = first(AG.filelist(dataset))
    (filepath=path, crs=crs, scale=scale, offset=offset, color=color, units=units)
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
