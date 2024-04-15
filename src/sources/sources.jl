# Source dispatch singletons
abstract type Source end

struct GRDsource <: Source end
struct GDALsource <: Source end

abstract type CDMsource <: Source end

struct GRIBsource <: CDMsource end
struct NCDsource <: CDMsource end

# Deprecations
const CDMfile = CDMsource
const NCDfile = NCDsource
const GRIBfile = GRIBsource
const GRDfile = GRDsource
const GDALfile = GDALsource

const SYMBOL2SOURCE = Dict(
    :gdal => GDALsource(),
    :grd => GRDsource(),
    :netcdf => NCDsource(), 
    :grib => GRIBsource(), 
)

const SOURCE2SYMBOL = Dict(map(reverse, collect(pairs(SYMBOL2SOURCE))))

# File extensions. GDAL is the catch-all for everything else
const SOURCE2EXT = Dict(
    GRDsource() => (".grd", ".gri"), 
    NCDsource() => (".nc", ".h5",), 
    GRIBsource() => (".grib",), 
)
const SOURCE2PACKAGENAME = Dict(
    GDALsource() => "ArchGDAL",
    NCDsource() => "NCDatasets",
    GRIBsource() => "GRIBDatasets",
)

const EXT2SOURCE = Dict(
    ".grd" => GRDsource(), 
    ".gri" => GRDsource(), 
    ".nc" => NCDsource(), 
    ".h5" => NCDsource(),
    ".grib" => GRIBsource(), 
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
_sourcetrait(filename::AbstractString, s::Source) = s
_sourcetrait(filename::AbstractString, s) = _sourcetrait(s)
_sourcetrait(filename::AbstractString, ::Union{Nothing,NoKW}) = _sourcetrait(filename)
_sourcetrait(filename::AbstractString) = get(EXT2SOURCE, splitext(filename)[2], GDALsource())
_sourcetrait(filenames::NamedTuple) = _sourcetrait(first(filenames))
_sourcetrait(filename, ext) = get(EXT2SOURCE, ext, GDALsource())
_sourcetrait(source::Source) = source
_sourcetrait(source::Type{<:Source}) = source()
function _sourcetrait(name::Symbol) 
    if haskey(SYMBOL2SOURCE, name)
        SYMBOL2SOURCE[name]
    else
        throw(ArgumentError("There is no source matching $name, try one from $(keys(SYMBOL2SOURCE))"))
    end
end

# Internal read method
function _open(f, filename::AbstractString; source=_sourcetrait(filename), kw...)
    _open(f, source, filename; kw...)
end
function _open(f, s::Source, filename::AbstractString; kw...)
    packagename = SOURCE2PACKAGENAME[s]
    throw(BackendException(packagename))
end
