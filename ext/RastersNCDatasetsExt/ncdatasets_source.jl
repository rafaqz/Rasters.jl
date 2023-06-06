const NCD = NCDatasets

"""
    Base.write(filename::AbstractString, ::Type{<:CDMsource}, A::AbstractRaster)

Write an NCDarray to a NetCDF file using NCDatasets.jl

Returns `filename`.
"""
function Base.write(filename::AbstractString, ::Type{<:CDMsource}, A::AbstractRaster; 
    append=false, force=false, verbose=true, kw...
)
    mode = if append
        isfile(filename) ? "a" : "c"
    else
        RA.check_can_write(filename, force)
        "c"
    end
    mode  = !isfile(filename) || !append ? "c" : "a";
    ds = NCD.Dataset(filename, mode; attrib=_attribdict(metadata(A)))
    try
        _ncdwritevar!(ds, A; kw...)
    finally
        close(ds)
    end
    return filename
end

# Stack ########################################################################

"""
    Base.write(filename::AbstractString, ::Type{NCDsource}, s::AbstractRasterStack; kw...)

Write an NCDstack to a single netcdf file, using NCDatasets.jl.

Currently `Metadata` is not handled for dimensions, and `Metadata` from other
[`AbstractRaster`](@ref) @types is ignored.

# Keywords

Keywords are passed to `NCDatasets.defVar`.

- `append`: If true, the variable of the current Raster will be appended to
    `filename`. Note that the variable of the current Raster should be not exist
    before. If not, you need to set `append = false`. Rasters.jl can not
    overwrite a previous existing variable.
- `fillvalue`: A value filled in the NetCDF file to indicate missing data. It
    will be stored in the `_FillValue` attribute.
- `chunksizes`: Vector integers setting the chunk size. The total size of a
    chunk must be less than 4 GiB.
- `deflatelevel`: Compression level: 0 (default) means no compression and 9
    means maximum compression. Each chunk will be compressed individually.
- `shuffle`: If true, the shuffle filter is activated which can improve the
    compression ratio.
- `checksum`: The checksum method can be `:fletcher32` or `:nochecksum`
    (checksumming is disabled, which is the default)
 - `typename` (string): The name of the NetCDF type required for vlen arrays
    (https://web.archive.org/save/https://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf-c/nc_005fdef_005fvlen.html)
"""
function Base.write(filename::AbstractString, ::Type{<:CDMsource}, s::AbstractRasterStack; append = false, kw...)
    mode  = !isfile(filename) || !append ? "c" : "a";
    ds = NCD.Dataset(filename, mode; attrib=_attribdict(metadata(s)))
    try
        map(key -> _ncdwritevar!(ds, s[key]), keys(s); kw...)
    finally
        close(ds)
    end
    return filename
end

RA.FileStack{NCDsource}(ds::AbstractDataset, filename::AbstractString; write=false, keys) = RA.FileStack(NCDsource, ds, filename; write, keys)

function RA.OpenStack(fs::RA.FileStack{NCDsource,K}) where K
    OpenStack{NCDsource,K}(NCD.Dataset(filename(fs)))
end
Base.close(os::RA.OpenStack{NCDsource}) = NCD.close(dataset(os))

function RA._open(f, ::Type{NCDsource}, filename::AbstractString; write=false, kw...)
    isfile(filename) || RA._isurl(filename) || RA._filenotfound_error(filename)
    mode = write ? "a" : "r"
    NCD.Dataset(filename, mode) do ds
        RA._open(f, NCDsource, ds; kw...)
    end
end

const _NCDVar = NCDatasets.CFVariable{Union{Missing, Float32}, 3, NCDatasets.Variable{Float32, 3, NCDatasets.NCDataset}, NCDatasets.Attributes{NCDatasets.NCDataset{Nothing}}, NamedTuple{(:fillvalue, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor), Tuple{Float32, Nothing, Nothing, Nothing, Nothing, Nothing}}}

# precompilation


# function _precompile(::Type{NCDsource})
#     ccall(:jl_generating_output, Cint, ()) == 1 || return nothing

#     precompile(Rasters.FileArray, (_NCDVar, String))
#     precompile(layerkeys, (NCDatasets.NCDataset{Nothing},))
#     precompile(dims, (_NCDVar,Symbol))
#     precompile(dims, (_NCDVar,Symbol,Nothing,Nothing))
#     precompile(dims, (_NCDVar,Symbol,Nothing,EPSG))
#     precompile(dims, (_NCDVar,Symbol,EPSG,EPSG))
#     precompile(_firstkey, (NCDatasets.NCDataset{Nothing},))
#     precompile(_ncddim, (NCDatasets.NCDataset{Nothing}, Symbol, Nothing, Nothing))
#     precompile(_ncddim, (NCDatasets.NCDataset{Nothing}, Symbol, Nothing, EPSG))
#     precompile(_ncddim, (NCDatasets.NCDataset{Nothing}, Symbol, EPSG, EPSG))
#     precompile(Raster, (NCDatasets.NCDataset{Nothing}, String, Nothing))
#     precompile(Raster, (NCDatasets.NCDataset{Nothing}, String, Symbol))
#     precompile(Raster, (_NCDVar, String, Symbol))

#     precompile(Raster, (String,))
# end

# _precompile(NCDsource)