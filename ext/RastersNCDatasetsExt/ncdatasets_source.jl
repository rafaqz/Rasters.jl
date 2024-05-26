const NCD = NCDatasets

const UNNAMED_NCD_FILE_KEY = "unnamed"

const NCDAllowedType = Union{Int8,UInt8,Int16,UInt16,Int32,UInt32,Int64,UInt64,Float32,Float64,Char,String}

function Base.write(filename::AbstractString, ::NCDsource, A::AbstractRaster;
    append=false,
    force=false,
    kw...
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
function Base.write(filename::AbstractString, ::NCDsource, s::AbstractRasterStack;
    append=false,
    force=false,
    missingval=nokw,
    kw...
)
    mode = if append
        isfile(filename) ? "a" : "c"
    else
        RA.check_can_write(filename, force)
        "c"
    end
    ds = NCD.Dataset(filename, mode; attrib=RA._attribdict(metadata(s)))
    try
        if missingval isa NamedTuple
            map(k -> _writevar!(ds, s[k]; missinval=missingval[k], kw...), keys(s))
        else
            map(k -> _writevar!(ds, s[k]; missingval, kw...), keys(s))
        end
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
function _writevar!(ds::AbstractDataset, A::AbstractRaster{T,N};
    verbose=true,
    missingval=nokw,
    chunks=nokw,
    chunksizes=RA._chunks_to_tuple(A, dims(A), chunks),
    scale=nokw,
    offset=nokw,
    eltype=Missings.nonmissingtype(T),
    kw...
) where {T,N}
    _def_dim_var!(ds, A)
    attrib = RA._attribdict(metadata(A))
    # Scale and offset
    scale = if isnokw(scale) || isnothing(scale)
        delete!(attrib, "scale_factor")
        nothing
    else
        attrib["scale_factor"] = scale
    end
    offset = if isnokw(offset) || isnothing(offset)
        delete!(attrib, "add_offset")
        nothing
    else
        attrib["add_offset"] = offset
    end
    maskingval1 = begin
        mv = maskingval isa NoKW ? Rasters.missingval(A) : maskingval
        mv === missingval ? nothing : mv
    end
    mod = _mod(missingval1, maskingval1, scale, offset)

    eltype <: NCDAllowedType || throw(ArgumentError("""
       Element type $eltyp cannot be written to NetCDF. Convert it to one of $(Base.uniontypes(NCDAllowedType)),
       usually by broadcasting the desired type constructor over the `Raster`, e.g. `newrast = Float32.(rast)`"))
       """
    ))

    # Set _FillValue
    if !isnothing(maskingval1) && Rasters.missingval(A) isa T
        attrib["_FillValue"] = missingval
    end
    verbose && !(maskingval isa Nothing) && @warn "`maskingval` $(maskingval) is not the same type as your data $T."

    key = if string(DD.name(A)) == ""
        UNNAMED_NCD_FILE_KEY
    else
        string(DD.name(A))
    end

    dimnames = lowercase.(string.(map(RA.name, dims(A))))
    var = _maybe_modify(NCD.defVar(ds, key, eltyp, dimnames; attrib=attrib, chunksizes, kw...), mod)
    # Write with a DiskArays.jl broadcast
    var .= A

    return nothing
end

_def_dim_var!(ds::AbstractDataset, A) = map(d -> _def_dim_var!(ds, d), dims(A))
function _def_dim_var!(ds::AbstractDataset, dim::Dimension)
    dimname = lowercase(string(DD.name(dim)))
    haskey(ds.dim, dimname) && return nothing
    NCD.defDim(ds, dimname, length(dim))
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
        boundskey = get(metadata(dim), :bounds, string(dimname, "_bnds"))
        push!(attrib, "bounds" => boundskey)
        NCD.defVar(ds, boundskey, bounds, ("bnds", dimname))
    end
    NCD.defVar(ds, dimname, Vector(index(dim)), (dimname,); attrib=attrib)
    return nothing
end

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
