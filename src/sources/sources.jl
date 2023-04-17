# Source dispatch singletons
abstract type Source end

abstract type CDMsource <: Source end
struct NCDsource <: CDMsource end
struct GRIBsource <: CDMsource end
struct GRDsource <: Source end
struct GDALsource <: Source end
struct SMAPsource <: Source end

# Deprecations
const NCDfile = NCDsource
const GRDfile = GRDsource
const GDALfile = GDALsource
const SMAPfile = SMAPsource

const SYMBOL2SOURCE = Dict(
    :gdal => GDALsource,
    :grd => GRDsource,
    :netcdf => NCDsource, 
    :smap => SMAPsource,
)

# File extensions. GDAL is the catch-all for everything else
const SOURCE2EXT = Dict(
    GRDsource => (".grd", ".gri"), 
    NCDsource => (".nc",), 
    SMAPsource => (".h5",)
)
const EXT2SOURCE = Dict(
    ".grd" => GRDsource, 
    ".gri" => GRDsource, 
    ".nc" => NCDsource, 
    ".grib" => GRIBsource, 
    ".h5" => SMAPsource
)

# Get the source backend for a file extension, falling back to GDALsource
_sourcetype(filename::AbstractString) = get(EXT2SOURCE, splitext(filename)[2], GDALsource)
_sourcetype(filenames::NamedTuple) = _sourcetype(first(filenames))
_sourcetype(filename, ext::Nothing) = _sourcetype(filename)
_sourcetype(filename, ext) = get(EXT2SOURCE, ext, GDALsource)
_sourcetype(source::Source) = typeof(source)
_sourcetype(source::Type{<:Source}) = source
function _sourcetype(name::Symbol) 
    if haskey(SYMBOL2SOURCE, name)
        SYMBOL2SOURCE[name]
    else
        throw(ArgumentError("There is no source matching $name, try one from $(keys(SYMBOL2SOURCE))"))
    end
end

# Internal read method
function _open(f, filename::AbstractString; source=_sourcetype(filename), kw...)
    _open(f, source, filename; kw...)
end
