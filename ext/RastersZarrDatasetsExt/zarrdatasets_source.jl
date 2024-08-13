function RA.OpenStack(fs::RA.FileStack{Zarrsource,K}) where K
    RA.OpenStack{Zarrsource,K}(ZD.ZarrDataset(RA.filename(fs)))
end

# In ZarrDatasets, the file is open for reading the values and closed afterwards. 
Base.close(os::RA.OpenStack{Zarrsource}) = nothing

function RA._open(f, ::Zarrsource, filename::AbstractString; write=false, kw...)
    ds = ZarrDatasets.ZarrDataset(filename)
    RA._open(f, Zarrsource(), ds; kw...)
end

RA._sourcetrait(::ZD.ZarrVariable) = Zarrsource()
RA._sourcetrait(::ZD.ZarrDataset) = Zarrsource()
