const NCDAllowedType = Union{Int8,UInt8,Int16,UInt16,Int32,UInt32,Int64,UInt64,Float32,Float64,Char,String}

function RA._check_allowed_type(::RA.NCDsource, eltyp)
    eltyp <: NCDAllowedType || throw(ArgumentError("""
    Element type $eltyp cannot be written to NetCDF. Convert it to one of $(Base.uniontypes(NCDAllowedType)),
    usually by broadcasting the desired type constructor over the `Raster`, e.g. `newrast = Float32.(rast)`"))
    """
    ))
end

function Base.write(filename::AbstractString, source::NCDsource, A::AbstractRaster;
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
        RA._writevar!(ds, source, A; kw...)
    finally
        close(ds)
    end
    return filename
end
function Base.write(filename::AbstractString, source::Source, s::AbstractRasterStack{K,T};
    append=false,
    force=false,
    missingval=nokw,
    maskingval=nokw,
    f=identity,
    kw...
) where {Source<:NCDsource,K,T}
    mode = if append
        isfile(filename) ? "a" : "c"
    else
        RA.check_can_write(filename, force)
        "c"
    end
    ds = NCD.Dataset(filename, mode; attrib=RA._attribdict(metadata(s)))

    maskingval = RA._stack_nt(s, isnokw(maskingval) ? Rasters.missingval(s) : maskingval)
    missingval = RA._stack_missingvals(s, isnokw(missingval) ? maskingval : missingval)
    try
        map(keys(s)) do k
            RA._writevar!(ds, source, s[k]; 
                missingval=missingval[k], 
                maskingval=maskingval[k], 
                kw...
            )
        end
        f(RA.OpenStack{Source,K,T}(ds))
    finally
        close(ds)
    end
    return filename
end

Base.close(os::RA.OpenStack{NCDsource}) = NCD.close(RA.dataset(os))

function RA.OpenStack(fs::RA.FileStack{NCDsource,K}) where K
    RA.OpenStack{NCDsource,K}(NCD.Dataset(RA.filename(fs)), fs.mods)
end

function RA._open(f, ::NCDsource, filename::AbstractString; write=false, kw...)
    isfile(filename) || RA._isurl(filename) || RA._filenotfound_error(filename)
    mode = write ? "a" : "r"
    NCD.Dataset(filename, mode) do ds
        RA._open(f, NCDsource(), ds; kw...)
    end
end

RA._sourcetrait(::NCD.Dataset) = NCDsource()
RA._sourcetrait(::NCD.Variable) = NCDsource()

@inline function RA.get_scale(metadata::Metadata{NCDsource}, scaled::Bool)
    scale = scaled ? get(metadata, "scale_factor", nothing) : nothing
    offset = scaled ? get(metadata, "add_offset", nothing) : nothing
    return scale, offset
end

RA.missingval(var::NCD.Variable, args...) = 
    RA.missingval(RA.Metadata{NCDsource}(CDM.attribs(var)))
RA.missingval(var::NCD.Variable, md::RA.Metadata{<:NCDsource}) = RA.missingval(md)

function RA.missingval(md::RA.Metadata{NCDsource})
    # TODO: handle multiple missing values
    fv = get(md, "_FillValue", nothing)
    mv = get(md, "missing_value", nothing)
    if isnothing(fv)
        if mv isa Vector
            length(mv) > 1 && @warn "'missing_value' $mv has multiple values. Currently we only uses the first."
            return first(mv)
        else
            return mv
        end
    else
        if !isnothing(mv) 
            fv == mv || @warn "Both '_FillValue' $fv and 'missing_value' $mv were found. Currently we only use the first."
        end
        return fv
    end
end

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
