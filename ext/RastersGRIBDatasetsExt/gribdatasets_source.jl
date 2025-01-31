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
# GRIB only has the parent object as a DiskArray
RA._open(f, ::GRIBsource, var::AbstractArray; mod=NoMod(), kw...) = 
    RA.cleanreturn(f(RA._maybe_modify(parent(var), mod)))

RA._sourcetrait(::GDS.Variable) = GRIBsource()
RA._sourcetrait(::GDS.GRIBDataset) = GRIBsource()

RA.missingval(var::GDS.Variable, ::RA.Metadata{GRIBsource}) = _missingval(var)
RA.missingval(var::GDS.Variable, args...) = _missingval(var)

function _missingval(var::GDS.Variable{T}) where T
    mv = GDS.missing_value(var)
    T1 = promote_type(typeof(mv), T)
    return T1(mv)
end
