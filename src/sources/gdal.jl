using ArchGDAL, GDAL, CoordinateTransformations

export GDALarray, GDALstack

# DimensionalData methods for ArchGDAL types ###############################

dims(dataset::ArchGDAL.Dataset, args...) = begin
    affine = geotransform_to_affine(ArchGDAL.getgeotransform(dataset))
    sze = ArchGDAL.width(dataset), ArchGDAL.height(dataset) 
    ax, ay = gdal_coord_convert(affine([0, 0]))
    bx, by = gdal_coord_convert(affine([sze...]))
    # println("ax, bx: ", (ax, bx))
    # println("ay, by: ", (ay, by))
    nbands = ArchGDAL.nraster(dataset)
    # TODO get an affine transform from the transformation
    lon = Lon((min(ax, bx), max(ax, bx))) 
    lat = Lat((min(ay, by), max(ay, by)); order=Order(Forward(), Reverse()))
    band = Band(1:nbands)
    formatdims((sze..., nbands), (lon, lat, band)) 
end

missingval(dataset::ArchGDAL.Dataset, args...) = 
    ArchGDAL.getnodatavalue(ArchGDAL.getband(dataset, 1))
metadata(dataset::ArchGDAL.Dataset, args...) = begin
    band = ArchGDAL.getband(dataset, 1)
    color = ArchGDAL.getname(ArchGDAL.getcolorinterp(band))
    scale = ArchGDAL.getscale(band)
    offset = ArchGDAL.getoffset(band)
    norvw = ArchGDAL.noverview(band)
    units = ArchGDAL.getunittype(band)
    crs = WellKnownText(string(ArchGDAL.getproj(dataset)))
    path = first(ArchGDAL.filelist(dataset))
    (filepath=path, crs=crs, scale=scale, offset=offset, color=color, units=units)
end

# Array ########################################################################

@GeoArrayMixin struct GDALarray{W,S,A} <: AbstractGeoArray{T,N,D}
    window::W
    size::S
end

GDALarray(path::AbstractString; kwargs...) = 
    gdalapply(dataset -> GDALarray(dataset; kwargs...), path)
GDALarray(dataset::ArchGDAL.Dataset; 
          dims=dims(dataset),
          refdims=(),
          metadata=metadata(dataset),
          missingval=missingval(dataset),
          name=Symbol(""),
          window=()) = begin
    path = first(ArchGDAL.filelist(dataset))
    if window == () 
        sze = gdalsize(dataset) 
    else
        sze = windowsize(window)
        dims, refdims = slicedims(dims, refdims, window)
    end
    T = ArchGDAL.getdatatype(ArchGDAL.getband(dataset, 1))
    N = length(sze)
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
    copy!(dst, gdalapply(ArchGDAL.read, source(src, key)))
Base.copy!(dst::AbstractGeoArray, src::GDALstack, key::Key) = 
copy!(parent(dst), gdalapply(ArchGDAL.read, source(src, key)))


# Utils ########################################################################

gdalapply(f, path::AbstractString) = 
    ArchGDAL.registerdrivers() do
        ArchGDAL.read(path) do dataset::ArchGDAL.Dataset
            f(dataset)
        end
    end

gdalread(dataset, I...) = ArchGDAL.read(dataset, reverse(I)...)
gdalread(dataset) = ArchGDAL.read(dataset)

gdalsize(dataset) = begin
    band = ArchGDAL.getband(dataset, 1)
    ArchGDAL.width(band), ArchGDAL.height(band), ArchGDAL.nraster(dataset)
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
get_affine_map(dataset::ArchGDAL.Dataset) = begin
    # ArchGDAL fails hard on datasets without
    # an affinemap. GDAL documents that on fail
    # a default affinemap should be returned.
    try
        global gt = ArchGDAL.getgeotransform(dataset)
    catch y
        @warn y.msg
        global gt = [0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
    end
    geotransform_to_affine(gt)
end

gdal_coord_convert(x) = x * 1e-5
