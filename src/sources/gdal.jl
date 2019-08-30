using ArchGDAL, GDAL, CoordinateTransformations

export GDALarray, GDALstack

# Interface methods for ArchGDAL types

@inline dims(ds::ArchGDAL.Dataset, args...) = begin
    affine = geotransform_to_affine(ArchGDAL.getgeotransform(ds))
    sze = ArchGDAL.width(ds), ArchGDAL.height(ds) 
    ax, ay = affine([0, 0])
    bx, by = affine([sze...])
    println((ax, bx))
    println((ay, by))
    nbands = ArchGDAL.nraster(ds)
    # TODO get an affine transform from the transformation
    lon = Lon <| (min(ax, bx), max(ax, bx)) 
    lat = Lat <| (min(ay, by), max(ay, by))
    band = Band(1:nbands)
    formatdims((sze..., nbands), (lon, lat, band)) 
end
@inline missingval(ds::ArchGDAL.Dataset, args...) = ArchGDAL.getnodatavalue(ArchGDAL.getband(ds, 1))
@inline metadata(ds::ArchGDAL.Dataset, args...) = begin
    band = ArchGDAL.getband(ds, 1)
    color = ArchGDAL.getname(ArchGDAL.getcolorinterp(band))
    scale = ArchGDAL.getscale(band)
    offset = ArchGDAL.getoffset(band)
    norvw = ArchGDAL.noverview(band)
    units = ArchGDAL.getunittype(band)
    crs = WellKnownText(string(ArchGDAL.importWKT(ArchGDAL.getproj(ds))))
    (filename=filename, crs=crs, scale=scale, offset=offset, color=color, units=units)
end


struct GDALarray{T,N,D,R,A,Me,Mi,S} <: AbstractGeoArray{T,N,D}
    path::A
    dims::D
    refdims::R
    metadata::Me
    missingval::Mi
    size::S
end
GDALarray(path::AbstractString; refdims=()) = 
    gdalrun(ds -> GDALarray(ds; refdims=refdims), path)
GDALarray(ds::ArchGDAL.Dataset; refdims=()) = begin
    dimz = dims(ds)
    meta = metadata(ds)
    mv = missingval(ds)
    path = first(ArchGDAL.filelist(ds))
    sze = gdalsize(ds)
    T = ArchGDAL.getdatatype(ArchGDAL.getband(ds, 1))
    N = length(sze)
    GDALarray{T,N,typeof.((dimz,refdims,path,meta,mv,sze))...
       }(path, dimz, refdims, meta, mv, sze)
end

Base.size(array::GDALarray) = array.size
Base.parent(array::GDALarray) = gdalrun(ds -> ArchGDAL.read(ds), array.path)
@inline Base.getindex(a::GDALarray, I::Vararg{<:Union{<:Integer,<:AbstractArray}}) = begin 
    data = gdalrun(ds -> ArchGDAL.read(ds, I...), a.path)
    GeoAarray(data, slicedims(a, I)..., metadata(a), missingval(s))
end

struct GDALstack{T} <: AbstractGeoStack{T} 
    data::T
end 

run(f, stack::GDALstack, path::AbstractString) = gdalrun(f, path::AbstractString)
data(stack::GDALstack, ds, key::Key, I::Vararg{Integer}) = ArchGDAL.read(ds, I...)
data(stack::GDALstack, ds, key::Key, I...) = GDALarray(ds)[I...] 
data(stack::GDALstack, ds, key::Key) = GDALarray(ds)

refdims(stack::GDALstack) = ()
metadata(stack::GDALstack) = nothing
metadata(stack::GDALstack, source, key::Key) = nothing
missingval(stack::GDALstack, source, key::Key) = NaN

crs(stack::GDALstack{<:NamedTuple}) = gdalrun(gdalcrs, first(parent(stack))[2])
crs(stack::GDALstack{<:AbstractString}) = gdalrun(gdalcrs, parent(stack))

Base.copy!(dst::AbstractArray, src::GDALstack, key::Key) = 
    copy!(dst, gdalrun(ArchGDAL.read, source(src, key)))


# gdal utils

gdalrun(f, path::AbstractString) = 
    ArchGDAL.registerdrivers(() -> ArchGDAL.read(f, path))

gdalsize(ds) = begin
    band = ArchGDAL.getband(ds, 1)
    ArchGDAL.width(band), ArchGDAL.height(band), ArchGDAL.nraster(ds)
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

# See https://lists.osgeo.org/pipermail/gdal-dev/2011-July/029449.html
# for an explanation of the geotransform format
geotransform_to_affine(gt) = AffineMap([gt[2] gt[3]; gt[5] gt[6]], [gt[1], gt[4]])
