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

GDALarray(filename::AbstractString; kwargs...) = begin
    isfile(filename) || error("file not found: $filename")
    gdalapply(dataset -> GDALarray(dataset; kwargs...), filename)
end
GDALarray(dataset::AG.Dataset;
          selectorcrs=nothing,
          dims=dims(dataset, selectorcrs),
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

Base.write(filename::AbstractString, ::Type{GDALarray},
           A::Union{<:GDALarray{T,2},<:GeoArray{T,2}}; kwargs...) where T = begin
    all(hasdim(A, (Lon, Lat))) || error("Array must have Lat and Lon dims to write to GTiff")
    A = permutedims(A, (Lon(), Lat()))
    nbands = 1
    indices = 1
    write_gdal(filename, A, nbands, indices)
end
Base.write(filename::AbstractString, ::Type{GDALarray},
           A::Union{<:GDALarray{T,3},<:GeoArray{T,3}}; kwargs...) where T = begin
    all(hasdim(A, (Lon, Lat))) || error("Array must have Lat and Lon dims to write to GeoTiff")
    hasdim(A, Band()) || error("Must have a `Band` dimension to write a 3-dimensional array to GeoTiff")
    A = permutedims(A, (Lon(), Lat(), Band()))
    nbands = size(A, Band())
    indices = Cint[1]
    write_gdal(filename, A, nbands, indices)
end

write_gdal(filename, A, nbands, indices; compress="DEFLATE", tiled="YES") = begin
    AG.create(filename;
        driver=AG.getdriver("GTiff"),
        width=size(A, 1),
        height=size(A, 2),
        nbands=nbands,
        dtype=eltype(A),
        options=["COMPRESS=$compress", "TILED=$tiled"],
       ) do dataset
        lon, lat = dims(A, (Lon(), Lat()))
        proj = convert(String, crs(grid(dims(A, Lat()))))
        AG.setproj!(dataset, proj)
        AG.setgeotransform!(dataset, build_geotransform(lat, lon))
        AG.write!(dataset, data(A), indices)
    end
    return filename
end

# TODO handle Forward/Reverse orders
build_geotransform(lat, lon) = begin
    gt = zeros(6)
    gt[GDAL_TOPLEFT_X] = first(lon)
    gt[GDAL_WE_RES] = step(lon)
    gt[GDAL_ROT1] = 0.0
    gt[GDAL_TOPLEFT_Y] = first(lat)
    gt[GDAL_ROT2] = 0.0
    gt[GDAL_NS_RES] = step(lat)
    return gt
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

dims(dataset::AG.Dataset, selectorcrs=nothing) = begin
    gt = try
        AG.getgeotransform(dataset)
    catch
        GDAL_EMPTY_TRANSFORM
    end

    latsize, lonsize = AG.height(dataset), AG.width(dataset)

    nbands = AG.nraster(dataset)
    band = Band(1:nbands, grid=CategoricalGrid())
    sourcecrs = crs(dataset)

    lonlat_metadata=GDALdimMetadata()

    # Output a BoundedGrid dims when the transformation is lat/lon alligned,
    # otherwise use TransformedGrid with an affine map.
    if isalligned(gt)
        lonstep = gt[GDAL_WE_RES]
        lonmin = gt[GDAL_TOPLEFT_X]
        lonmax = lonmin + lonstep * (lonsize - 1)
        lonrange = LinRange(lonmin, lonmax, lonsize)

        latstep = gt[GDAL_NS_RES]
        latmax = gt[GDAL_TOPLEFT_Y]
        latmin = latmax + latstep * (latsize - 1)
        latrange = LinRange(latmax, latmin, latsize)

        areaorpoint = gdalmetadata(dataset, "AREA_OR_POINT")
        # Spatial data defaults to area/inteval?
        if areaorpoint == "Point"
            sampling = PointSampling()
        else areaorpoint 
            # GeoTiff uses the "pixelCorner" convention
            sampling = IntervalSampling(Start())
        end

        latgrid = ProjectedGrid(
            # Latitude is in reverse to how plot it.
            order=Ordered(Reverse(), Reverse(), Forward()), 
            sampling=sampling,
            # Use the range step as is will be different to latstep due to float error
            span=RegularSpan(step(latrange)), 
            crs=sourcecrs, 
            selectorcrs=selectorcrs
        )
        longrid = ProjectedGrid(
            span=RegularSpan(step(lonrange)), 
            sampling=sampling,
            crs=sourcecrs, 
            selectorcrs=selectorcrs
        )

        lon = Lon(lonrange; grid=longrid, metadata=lonlat_metadata)
        lat = Lat(latrange; grid=latgrid, metadata=lonlat_metadata)

        formatdims(map(Base.OneTo, (lonsize, latsize, nbands)), (lon, lat, band))
    else
        error("Rotated/transformed grids not handled currently")
        # affinemap = geotransform_to_affine(geotransform)
        # x = X(affinemap; grid=TransformedGrid(dims=Lon()))
        # y = Y(affinemap; grid=TransformedGrid(dims=Lat()))

        # formatdims((lonsize, latsize, nbands), (x, y, band))
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

gdalmetadata(dataset, key) = begin
    meta = AG.metadata(dataset)
    regex = Regex("$key=(.*)")
    i = findfirst(f -> occursin(regex, f), meta)
    if i isa Nothing
        nothing
    else
        match(regex, meta[i])[1]
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

isalligned(geotransform) = 
    geotransform[GDAL_ROT1] == 0 && geotransform[GDAL_ROT2] == 0

geotransform_to_affine(gt) = begin
    AffineMap([gt[GDAL_WE_RES] gt[GDAL_ROT1]; gt[GDAL_ROT2] gt[GDAL_NS_RES]], 
              [gt[GDAL_TOPLEFT_X], gt[GDAL_TOPLEFT_Y]])
end


