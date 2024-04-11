const NCD = NCDatasets

const UNNAMED_NCD_FILE_KEY = "unnamed"

const NCDAllowedType = Union{Int8,UInt8,Int16,UInt16,Int32,UInt32,Int64,UInt64,Float32,Float64,Char,String}

"""
    Base.write(filename::AbstractString, ::NCDsource, A::AbstractRaster)

Write an NCDarray to a NetCDF file using NCDatasets.jl

Returns `filename`.
"""
function Base.write(filename::AbstractString, ::NCDsource, A::AbstractRaster; 
    append=false, force=false, verbose=true, kw...
)
    mode = if append
        isfile(filename) ? "a" : "c"
    else
        RA.check_can_write(filename, force)
        "c"
    end
    mode  = !isfile(filename) || !append ? "c" : "a";
    ds = NCD.Dataset(filename, mode; attrib=RA._attribdict(metadata(A)))
    try
        _writevar!(ds, A; kw...)
    finally
        close(ds)
    end
    return filename
end

# Stack ########################################################################

"""
    Base.write(filename::AbstractString, ::NCDsource, s::AbstractRasterStack; kw...)

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
function Base.write(filename::AbstractString, ::NCDsource, s::AbstractRasterStack; append = false, kw...)
    mode  = !isfile(filename) || !append ? "c" : "a";
    ds = NCD.Dataset(filename, mode; attrib=RA._attribdict(metadata(s)))
    try
        map(key -> _writevar!(ds, s[key]), keys(s); kw...)
    finally
        close(ds)
    end
    return filename
end

function RA.OpenStack(fs::RA.FileStack{NCDsource,K}) where K
    RA.OpenStack{NCDsource,K}(NCD.Dataset(RA.filename(fs)))
end
Base.close(os::RA.OpenStack{NCDsource}) = NCD.close(RA.dataset(os))

function RA._open(f, ::NCDsource, filename::AbstractString; write=false, kw...)
    isfile(filename) || RA._isurl(filename) || RA._filenotfound_error(filename)
    mode = write ? "a" : "r"
    NCD.Dataset(filename, mode) do ds
        RA._open(f, NCDsource(), ds; kw...)
    end
end

# Add a var array to a dataset before writing it.
function _writevar!(ds::AbstractDataset, A::AbstractRaster{T,N}; kw...) where {T,N}
    _def_dim_var!(ds, A)
    attrib = RA._attribdict(metadata(A))
    # Set _FillValue
    eltyp = Missings.nonmissingtype(T)
    eltyp <: NCDAllowedType || throw(ArgumentError("$eltyp cannot be written to NetCDF, convert to one of $(Base.uniontypes(NCDAllowedType))"))
    if ismissing(missingval(A))
        fillval = if haskey(attrib, "_FillValue") && attrib["_FillValue"] isa eltyp
            attrib["_FillValue"]
        else
            NCD.fillvalue(eltyp)
        end
        attrib["_FillValue"] = fillval
        A = replace_missing(A, fillval)
    elseif missingval(A) isa T
        attrib["_FillValue"] = missingval(A)
    else
        missingval(A) isa Nothing || @warn "`missingval` $(missingval(A)) is not the same type as your data $T."
    end

    key = if string(DD.name(A)) == ""
        UNNAMED_NCD_FILE_KEY
    else
        string(DD.name(A))
    end

    dimnames = lowercase.(string.(map(RA.name, dims(A))))
    var = NCD.defVar(ds, key, eltyp, dimnames; attrib=attrib, kw...) |> RA.CFDiskArray

    # Write with a DiskArays.jl broadcast
    var .= A

    return nothing
end

_def_dim_var!(ds::AbstractDataset, A) = map(d -> _def_dim_var!(ds, d), dims(A))
function _def_dim_var!(ds::AbstractDataset, dim::Dimension)
    dimkey = lowercase(string(DD.name(dim)))
    haskey(ds.dim, dimkey) && return nothing
    NCD.defDim(ds, dimkey, length(dim))
    lookup(dim) isa NoLookup && return nothing

    # Shift index before conversion to Mapped
    dim = RA._cdmshiftlocus(dim)
    if dim isa Y || dim isa X
        dim = convertlookup(Mapped, dim)
    end
    # Attributes
    attrib = RA._attribdict(metadata(dim))
    RA._cdm_set_axis_attrib!(attrib, dim)
    # Bounds variables
    if sampling(dim) isa Intervals
        bounds = Dimensions.dim2boundsmatrix(dim)
        boundskey = get(metadata(dim), :bounds, string(dimkey, "_bnds"))
        push!(attrib, "bounds" => boundskey)
        NCD.defVar(ds, boundskey, bounds, ("bnds", dimkey))
    end
    NCD.defVar(ds, dimkey, Vector(index(dim)), (dimkey,); attrib=attrib)
    return nothing
end

# Hack to get the inner DiskArrays chunks as they are not exposed at the top level
RA._get_eachchunk(var::NCD.Variable) = DiskArrays.eachchunk(var)
RA._get_haschunks(var::NCD.Variable) = DiskArrays.haschunks(var)

RA._sourcetrait(::NCD.Dataset) = NCDsource()
RA._sourcetrait(::NCD.Variable) = NCDsource()

# precompilation

# const _NCDVar = NCDatasets.CFVariable{Union{Missing, Float32}, 3, NCDatasets.Variable{Float32, 3, NCDatasets.NCDataset}, NCDatasets.Attributes{NCDatasets.NCDataset{Nothing}}, NamedTuple{(:fillvalue, :scale_factor, :add_offset, :calendar, :time_origin, :time_factor), Tuple{Float32, Nothing, Nothing, Nothing, Nothing, Nothing}}}

# function _precompile(::Type{NCDsource})
#     ccall(:jl_generating_output, Cint, ()) == 1 || return nothing

#     precompile(Rasters.FileArray, (_NCDVar, String))
#     precompile(layerkeys, (NCDatasets.NCDataset{Nothing},))
#     precompile(dims, (_NCDVar,Symbol))
#     precompile(dims, (_NCDVar,Symbol,Nothing,Nothing))
#     precompile(dims, (_NCDVar,Symbol,Nothing,EPSG))
#     precompile(dims, (_NCDVar,Symbol,EPSG,EPSG))
#     precompile(_firstkey, (NCDatasets.NCDataset{Nothing},))
#     precompile(_cdmdim, (NCDatasets.NCDataset{Nothing}, Symbol, Nothing, Nothing))
#     precompile(_cdmdim, (NCDatasets.NCDataset{Nothing}, Symbol, Nothing, EPSG))
#     precompile(_cdmdim, (NCDatasets.NCDataset{Nothing}, Symbol, EPSG, EPSG))
#     precompile(Raster, (NCDatasets.NCDataset{Nothing}, String, Nothing))
#     precompile(Raster, (NCDatasets.NCDataset{Nothing}, String, Symbol))
#     precompile(Raster, (_NCDVar, String, Symbol))

#     precompile(Raster, (String,))
# end

# _precompile(NCDsource)
