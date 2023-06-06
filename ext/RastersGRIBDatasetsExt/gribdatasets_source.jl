const GDS = GRIBDatasets

RA.FileStack{GRIBsource}(ds::AbstractDataset, filename::AbstractString; write=false, keys) = RA.FileStack(GRIBsource, ds, filename; write, keys)

function RA.OpenStack(fs::RA.FileStack{GRIBsource,K}) where K
    RA.OpenStack{GRIBsource,K}(GDS.GRIBDataset(RA.filename(fs)))
end

# In GRIBDatasets, the file is open for reading the values and closed afterwards. 
Base.close(os::RA.OpenStack{GRIBsource}) = nothing

function RA._open(f, ::Type{GRIBsource}, filename::AbstractString; write=false, kw...)
    isfile(filename) || _filenotfound_error(filename)
    ds = GRIBDatasets.GRIBDataset(filename)
    RA._open(f, GRIBsource, ds; kw...)
end
