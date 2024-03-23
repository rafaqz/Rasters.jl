# Source dispatch singletons
abstract type Source end

abstract type CDMsource <: Source end
struct NCDsource <: CDMsource end
struct GRIBsource <: CDMsource end
struct GRDsource <: Source end
struct GDALsource <: Source end
struct SMAPsource <: Source end

# Deprecations
const CDMfile = CDMsource
const NCDfile = NCDsource
const GRIBfile = GRIBsource
const GRDfile = GRDsource
const GDALfile = GDALsource
const SMAPfile = SMAPsource

const SYMBOL2SOURCE = Dict(
    :gdal => GDALsource(),
    :grd => GRDsource(),
    :netcdf => NCDsource(), 
    :grib => GRIBsource(), 
    :smap => SMAPsource(),
)

const SOURCE2SYMBOL = Dict(map(reverse, collect(pairs(SYMBOL2SOURCE))))

# File extensions. GDAL is the catch-all for everything else
const SOURCE2EXT = Dict(
    GRDsource() => (".grd", ".gri"), 
    NCDsource() => (".nc",), 
    GRIBsource() => (".grib",), 
    SMAPsource() => (".h5",),
)
const SOURCE2PACAKGENAME = Dict(
    GDALsource() => "ArchGDAL",
    NCDsource() => "NCDatasets",
    GRIBsource() => "GRIBDatasets",
    SMAPsource() => "HDF5",
)

const EXT2SOURCE = Dict(
    ".grd" => GRDsource(), 
    ".gri" => GRDsource(), 
    ".nc" => NCDsource(), 
    ".grib" => GRIBsource(), 
    ".h5" => SMAPsource(),
)

# exception to be raised when backend extension is not satisfied
struct BackendException <: Exception
    backend
end

# error message to show when backend is not loaded
function Base.showerror(io::IO, e::BackendException)
    print(io, "`Rasters.jl` requires backends to be loaded externally as of version 0.8. Run `import $(e.backend)` to fix this error.")
end

# Get the source backend for a file extension, falling back to GDALsource
_sourcetype(filename::AbstractString) = get(EXT2SOURCE, splitext(filename)[2], GDALsource())
_sourcetype(filenames::NamedTuple) = _sourcetype(first(filenames))
_sourcetype(filename, ext::Nothing) = _sourcetype(filename)
_sourcetype(filename, ext) = get(EXT2SOURCE, ext, GDALsource())
_sourcetype(source::Source) = source
_sourcetype(source::Type{<:Source}) = source()
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
function _open(f, s::Source, filename::AbstractString; kw...)
    packagename = SOURCE2PACAKGENAME[s]
    throw(BackendException(packagename))
end
