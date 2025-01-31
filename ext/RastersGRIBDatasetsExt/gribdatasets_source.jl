RA.sourceconstructor(::Type{GRIBsource}) = GDS.GRIBDataset

RA.checkmode(::GRIBsource, filename, append::Bool, force::Bool) =
    throw(ArgumentError("GRIBDatasets.jl does not support writing"))

# In GRIBDatasets, the file is open for reading the values and closed afterwards. 
Base.close(os::RA.OpenStack{GRIBsource}) = nothing

function RA._open(f, ::GRIBsource, filename::AbstractString; write=false, kw...)
    isfile(filename) || _filenotfound_error(filename)
    ds = GRIBDatasets.GRIBDataset(filename)
    RA._open(f, GRIBsource(), ds; kw...)
end

@inline function RA.get_scale(metadata::Metadata{NCDsource}, scaled::Bool)
    scale = scaled ? get(metadata, "scale_factor", nothing) : nothing
    offset = scaled ? get(metadata, "add_offset", nothing) : nothing
    return scale, offset
end

RA._sourcetrait(::GDS.Variable) = GRIBsource()
RA._sourcetrait(::GDS.GRIBDataset) = GRIBsource()

function RA.missingval(var::GDS.Variable{T}, args...) where T
    mv = GDS.missing_value(var)
    T1 = promote_type(typeof(mv), T)
    return T1(mv)
end
