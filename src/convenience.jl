# File extensions. GDAL is the catch-all for everything else
const EXT = Dict(
    GRDfile => (".grd", ".gri"), 
    NCDfile => (".nc",), 
    SMAPfile => (".h5",)
)
const REV_EXT = Dict(
    ".grd" => GRDfile, 
    ".gri" => GRDfile, 
    ".nc" => NCDfile, 
    ".h5" => SMAPfile,
    ".asc" => ASCIIfile
)

# Get the source backend for a file extension, falling back to GDALfile
_sourcetype(filename::AbstractString) = get(REV_EXT, splitext(filename)[2], GDALfile)
_sourcetype(filenames::NamedTuple) = _sourcetype(first(filenames))

# Internal read method
function _open(f, filename::AbstractString; kw...)
    _open(f, _sourcetype(filename), filename; kw...)
end
