const GDS = GRIBDatasets

function RA.OpenStack(fs::RA.FileStack{GRIBsource,K}) where K
    RA.OpenStack{GRIBsource,K}(GDS.GRIBDataset(RA.filename(fs)))
end

# In GRIBDatasets, the file is open for reading the values and closed afterwards. 
Base.close(os::RA.OpenStack{GRIBsource}) = nothing

function RA._open(f, ::GRIBsource, filename::AbstractString; write=false, kw...)
    isfile(filename) || _filenotfound_error(filename)
    ds = GRIBDatasets.GRIBDataset(filename)
    RA._open(f, GRIBsource(), ds; kw...)
end

# Hack to get the inner DiskArrays chunks as they are not exposed at the top level
RA._get_eachchunk(var::GDS.Variable) = DiskArrays.eachchunk(var.values)
RA._get_haschunks(var::GDS.Variable) = DiskArrays.haschunks(var.values)

RA._sourcetrait(::GDS.Variable) = GRIBsource()
RA._sourcetrait(::GDS.GRIBDataset) = GRIBsource()
