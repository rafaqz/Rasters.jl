function RA.OpenStack(fs::RA.FileStack{GRIBsource,K}) where K
    RA.OpenStack{GRIBsource,K}(GDS.GRIBDataset(RA.filename(fs)), fs.mods)
end

# In GRIBDatasets, the file is open for reading the values and closed afterwards. 
Base.close(os::RA.OpenStack{GRIBsource}) = nothing

function RA._open(f, ::GRIBsource, filename::AbstractString; write=false, kw...)
    isfile(filename) || _filenotfound_error(filename)
    ds = GRIBDatasets.GRIBDataset(filename)
    RA._open(f, GRIBsource(), ds; kw...)
end

RA._sourcetrait(::GDS.Variable) = GRIBsource()
RA._sourcetrait(::GDS.GRIBDataset) = GRIBsource()

function RA.missingval(var::GDS.Variable{T}, args...) where T
    mv = GDS.missing_value(var)
    T1 = promote_type(typeof(mv), T)
    return T1(mv)
end
