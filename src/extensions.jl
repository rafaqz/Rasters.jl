# extensions
function throw_extention_error(package::String, extention::Symbol)
    @static if isdefined(Base, :get_extension) # julia > 1.9
    if isnothing(Base.get_extension(Rasters, extention))
        throw(BackendException(package))
    else
        throw(MethodError(resample, args))
    end
    else
        throw(BackendException(package))
    end
end


# stubs that need ArchGDAL
resample(args...; kw...) = throw_extention_error("ArchGDAL", :RastersArchGDALExt)
warp(args...; kw...) = throw_extention_error("ArchGDAL", :RastersArchGDALExt)
cellsize(args...; kw...) = throw_extention_error("ArchGDAL", :RastersArchGDALExt)

# Other shared stubs
function layerkeys end
function smapseries end
function maybe_correct_to_write end
function dims2geotransform end
function affine2geotransform end
function geotransform2affine end

# Shared between ArchGDAL and CoordinateTransformations extenstions
const GDAL_EMPTY_TRANSFORM = [0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
const GDAL_TOPLEFT_X = 1
const GDAL_WE_RES = 2
const GDAL_ROT1 = 3
const GDAL_TOPLEFT_Y = 4
const GDAL_ROT2 = 5
const GDAL_NS_RES = 6

# The rest of the definition is in CoordinateTransformations
struct AffineProjected{T,F,A<:AbstractVector{T},M,C,MC,P,D} <: LA.Unaligned{T,1}
    affinemap::F
    data::A
    metadata::M
    crs::C
    mappedcrs::MC
    paired_lookup::P
    dim::D
end
