RA.sourcetrait(::GDS.Variable) = GRIBsource()
RA.sourcetrait(::GDS.GRIBDataset) = GRIBsource()

RA.sourceconstructor(::Type{GRIBsource}) = _gribdataset
# GRIB doesn't accept the mode keyword so hack around it
_gribdataset(filename, mode="") = GDS.GRIBDataset(filename)

# GribDatasets.jl is essentially broken to 
# use directly, so we get the internal values
# RA._open(f, ::GRIBsource, var::AbstractArray; mod=NoMod(), kw...) = 
    # RA.cleanreturn(f(RA._maybe_modify(var.values, mod)))

RA.checkwritemode(::GRIBsource, filename, append::Bool, force::Bool) =
    throw(ArgumentError("GRIBDatasets.jl does not support writing"))
RA.checkfilename(::GRIBsource, filename) =
    isfile(filename) || _filenotfound_error(filename)

# In GRIBDatasets, the file is open for reading the values and closed afterwards. 
Base.close(os::RA.OpenStack{GRIBsource}) = nothing

RA.missingval(var::GDS.Variable, ::RA.Metadata{<:RA.CDMsource}) = _missingval(var)
RA.missingval(var::GDS.Variable, args...) = _missingval(var)

function _missingval(var::GDS.Variable)
    mv = GDS.missing_value(var)
    RA._fix_missingval(var, mv)
end