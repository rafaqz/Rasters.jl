# Source dispatch singletons
"""
    Source
    
Abstract type for all sources.  This is used to dispatch on the source
backend to use for reading a file.  The source is determined by the file
extension, or by the source keyword argument.
"""
abstract type Source end

abstract type CDMsource <: Source end

struct GRIBsource <: CDMsource end
struct NCDsource <: CDMsource end
struct Zarrsource <: CDMsource end

struct GRDsource <: Source end
struct GDALsource <: Source end

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
    :zarr => Zarrsource(),
)

const SOURCE2SYMBOL = Dict(map(reverse, collect(pairs(SYMBOL2SOURCE))))

# File extensions. GDAL is the catch-all for everything else
const SOURCE2EXT = Dict(
    GRDsource() => (".grd", ".gri"), 
    NCDsource() => (".nc", ".nc4", ".h5",), 
    GRIBsource() => (".grib",), 
    Zarrsource() => (".zarr", ".zarr/"),
)
const SOURCE2PACKAGENAME = Dict(
    GDALsource() => "ArchGDAL",
    NCDsource() => "NCDatasets",
    GRIBsource() => "GRIBDatasets",
    Zarrsource() => "ZarrDatasets",
)

const EXT2SOURCE = Dict(
    ".grd" => GRDsource(), 
    ".gri" => GRDsource(), 
    ".nc" => NCDsource(), 
    ".nc4" => NCDsource(), 
    ".h5" => NCDsource(),
    ".grib" => GRIBsource(), 
    ".zarr" => Zarrsource(),
)

# exception to be raised when backend extension is not satisfied
struct BackendException <: Exception
    backend
end
BackendException(s::Source) = BackendException(SOURCE2PACKAGENAME[s])

# error message to show when backend is not loaded
function Base.showerror(io::IO, e::BackendException)
    printstyled(io, "Rasters.jl"; underline = true)
    printstyled(io, " requires backends to be loaded manually.  Run ")
    printstyled(io, "`import $(e.backend)`"; bold = true)
    print(io, " to fix this error.")
end

# Dataset constructor from `Source`
sourceconstructor(s::Source) = throw(BackendException(s))
# Function to check filename
checkfilename(s::Source, filename) = throw(BackendException(s))

# Get the source backend for a file extension, falling back to GDALsource
sourcetrait(x) = 
    throw(ArgumentError("`sourcetrait` is not defined for $(typeof(x))"))
sourcetrait(filename::AbstractString, s::Symbol) = sourcetrait(s)
sourcetrait(filename::AbstractString, s::Source) = s
sourcetrait(filename::AbstractString, ::Type{S}) where S<:Source = S()
sourcetrait(filename::AbstractString, ::Union{Nothing,NoKW}) = sourcetrait(filename)
sourcetrait(filename::AbstractString, ext::AbstractString) = get(EXT2SOURCE, ext, GDALsource())
function sourcetrait(filename::AbstractString)
    isempty(filename) && throw(ArgumentError("Filename cannot be empty"))
    default = GDALsource()
    stem, ext = splitext(filename)
    str = if ext == ""  
        # Handle e.g. "x.zarr/" directories
        if isdirpath(stem) 
            return sourcetrait(dirname(stem))
        else
            stem
        end
    else
        ext
    end
    return get(EXT2SOURCE, str, default)
end
sourcetrait(filenames::NamedTuple) = sourcetrait(first(filenames))
sourcetrait(source::Source) = source
sourcetrait(source::Type{<:Source}) = source()
function sourcetrait(name::Symbol) 
    if haskey(SYMBOL2SOURCE, name)
        SYMBOL2SOURCE[name]
    else
        throw(ArgumentError("There is no source matching $name, try one from $(keys(SYMBOL2SOURCE))"))
    end
end

# Internal read method
function _open(f, filename::AbstractString; source=sourcetrait(filename), kw...)
    _open(f, source, filename; kw...)
end

# Open a backend file and return the dataset handle (which may itself be a
# stack-like container or an array, depending on the backend). Backends
# define this. Used by both the non-closure `Base.open(::AbstractRaster)`
# path and the closure-form `_open` (via try/finally).
_open_dataset(s::Source, filename::AbstractString; kw...) = throw(BackendException(s))

# Open a single array from an already-open dataset. Returns the value to
# use as the parent of an opened `Raster` — typically the backend array
# wrapped in `_maybe_modify`. Backends define this.
_open_array(s::Source, ds; kw...) = throw(BackendException(s))

# Extract the dataset handle backing an opened value (CFVariable,
# RasterDataset, wrapper, …) or `nothing` if there is no resource to close
# (in-memory arrays, GRD mmaps). Backends extend with methods on their
# array/dataset types.
_dataset(::Any) = nothing
# Generic peel: a DiskArray wrapper with a distinct parent forwards to it.
# Concrete leaves (CFVariable, RasterDataset, RasterDiskArray over an mmap)
# add their own methods; this handles intermediate wrappers like
# `DiskArrays.SubDiskArray` and `DiskArrays.BroadcastDiskArray` without
# having to enumerate each one.
function _dataset(A::DiskArrays.AbstractDiskArray)
    p = Base.parent(A)
    return p === A ? nothing : _dataset(p)
end

# Close a dataset handle. Default is a no-op for stateless backends (Zarr,
# GRIB, GRD). Backends with real handles override on their dataset type.
_close_dataset(::Any) = nothing

# Closure-form `_open` built on top of `_open_dataset` / `_close_dataset`.
# Backends with bespoke lifetimes (e.g. GRD mmap) can still override.
function _open(f, source::Source, filename::AbstractString; kw...)
    ds = _open_dataset(source, filename; kw...)
    try
        return _open(f, source, ds; kw...)
    finally
        _close_dataset(ds)
    end
end
