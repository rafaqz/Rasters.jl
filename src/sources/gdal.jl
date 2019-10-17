using ArchGDAL, GDAL, CoordinateTransformations

export GDALarray, GDALstack

# DimensionalData methods for ArchGDAL types ###############################

dims(ds::ArchGDAL.Dataset, args...) = begin
    affine = geotransform_to_affine(ArchGDAL.getgeotransform(ds))
    sze = ArchGDAL.width(ds), ArchGDAL.height(ds) 
    ax, ay = gdal_coord_convert(affine([0, 0]))
    bx, by = gdal_coord_convert(affine([sze...]))
    # println("ax, bx: ", (ax, bx))
    # println("ay, by: ", (ay, by))
    nbands = ArchGDAL.nraster(ds)
    # TODO get an affine transform from the transformation
    lon = Lon((min(ax, bx), max(ax, bx))) 
    lat = Lat((min(ay, by), max(ay, by)); order=Order(Forward(), Reverse()))
    band = Band(1:nbands)
    formatdims((sze..., nbands), (lon, lat, band)) 
end

missingval(ds::ArchGDAL.Dataset, args...) = ArchGDAL.getnodatavalue(ArchGDAL.getband(ds, 1))
metadata(ds::ArchGDAL.Dataset, args...) = begin
    band = ArchGDAL.getband(ds, 1)
    color = ArchGDAL.getname(ArchGDAL.getcolorinterp(band))
    scale = ArchGDAL.getscale(band)
    offset = ArchGDAL.getoffset(band)
    norvw = ArchGDAL.noverview(band)
    units = ArchGDAL.getunittype(band)
    crs = WellKnownText(string(ArchGDAL.getproj(ds)))
    path = first(ArchGDAL.filelist(ds))
    (filepath=path, crs=crs, scale=scale, offset=offset, color=color, units=units)
end

# Array ########################################################################

@GeoArrayMixin struct GDALarray{W,S,A} <: AbstractGeoArray{T,N,D}
    window::W
    size::S
end

GDALarray(path::AbstractString; kwargs...) = gdalapply(ds -> GDALarray(ds; kwargs...), path)
GDALarray(ds::ArchGDAL.Dataset; 
          dims=dims(ds),
          refdims=(),
          metadata=metadata(ds),
          missingval=missingval(ds),
          name=Symbol(""),
          window=()) = begin
    path = first(ArchGDAL.filelist(ds))
    if window == () 
        sze = gdalsize(ds) 
    else
        sze = windowsize(window)
        dims, refdims = slicedims(dims, refdims, window)
    end
    T = ArchGDAL.getdatatype(ArchGDAL.getband(ds, 1))
    N = length(sze)
    GDALarray{T,N,typeof.((dims,refdims,metadata,missingval,name,window,sze,path))...
       }(path, dims, refdims, metadata, missingval, name, window, sze)
end

Base.size(a::GDALarray) = a.size
Base.parent(a::GDALarray) = gdalapply(ds -> gdalread(ds, windoworempty(a)...), a.data)
Base.getindex(a::GDALarray, I::Vararg{<:Union{<:Integer,<:AbstractArray}}) = begin
    I = applywindow(a, I)
    rebuildsliced(a, gdalapply(ds -> gdalread(ds, I...), a.data), I)
end

# Stack ########################################################################

@GeoStackMixin struct GDALstack{} <: AbstractGeoStack{T} end 

GDALstack(data::NamedTuple; 
          dims=gdalapply(dims, first(values(data))), 
          refdims=(), window=(), 
          metadata=gdalapply(metadata, first(values(data)))) =
    GDALstack(data, dims, refdims, window, metadata)

safeapply(f, ::GDALstack, path::AbstractString) = gdalapply(f, path)
data(stack::GDALstack, ds, key::Key, I...) = GDALarray(ds; window=I)
data(stack::GDALstack, ds, key::Key) = GDALarray(ds)

Base.copy!(dst::AbstractArray, src::GDALstack, key::Key) = 
    copy!(dst, gdalapply(ArchGDAL.read, source(src, key)))


# Utils ########################################################################

gdalapply(f, path::AbstractString) = 
    ArchGDAL.registerdrivers() do
        ArchGDAL.read(path) do ds::ArchGDAL.Dataset
            f(ds)
        end
    end

gdalread(ds, I...) = ArchGDAL.read(ds, reverse(I)...)
gdalread(ds) = ArchGDAL.read(ds)

gdalsize(ds) = begin
    band = ArchGDAL.getband(ds, 1)
    ArchGDAL.width(band), ArchGDAL.height(band), ArchGDAL.nraster(ds)
end

# See https://lists.osgeo.org/pipermail/gdal-dev/2011-July/029449.html
# for an explanation of the geotransform format
geotransform_to_affine(gt) = begin
    # println(gt)
    AffineMap([gt[2] gt[3]; gt[5] gt[6]], [gt[1], gt[4]])
end

# This is copied from evetion/GeoArrays.jl
# It's not yet clear how to integrate these packages so
# so this is really a placeholder for that happening.
get_affine_map(ds::ArchGDAL.Dataset) = begin
    # ArchGDAL fails hard on datasets without
    # an affinemap. GDAL documents that on fail
    # a default affinemap should be returned.
    try
        global gt = ArchGDAL.getgeotransform(ds)
    catch y
        @warn y.msg
        global gt = [0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
    end
    geotransform_to_affine(gt)
end

gdal_coord_convert(x) = x * 1e-5
