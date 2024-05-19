function RA.OpenStack(fs::RA.FileStack{Zarrsource,K}) where K
    RA.OpenStack{Zarrsource,K}(ZD.ZarrDataset(RA.filename(fs)))
end

# In ZarrDatasets, the file is open for reading the values and closed afterwards. 
Base.close(os::RA.OpenStack{Zarrsource}) = nothing

function RA._open(f, ::Zarrsource, filename::AbstractString; write=false, kw...)
    #isfile(filename) || RA._filenotfound_error(filename)
    ds = ZarrDatasets.ZarrDataset(filename)
    RA._open(f, Zarrsource(), ds; kw...)
end

function Base.write(filename::AbstractString, ::Zarrsource, A::AbstractRaster;
    append=false,
    force=false,
    kw...
)
    writeable = RA.check_can_write(filename, force)
    mode="c"
    ds = ZarrDataset(filename, mode; attrib=RA._attribdict(metadata(A)))
    try
        RA._writevar!(ds, A; kw...)
    finally
        close(ds)
    end
    return filename
end


# Hack to get the inner DiskArrays chunks as they are not exposed at the top level
RA._get_eachchunk(var::ZD.ZarrVariable) = DiskArrays.eachchunk(var.zarray)
RA._get_haschunks(var::ZD.ZarrVariable) = DiskArrays.haschunks(var.zarray)

RA._sourcetrait(::ZD.ZarrVariable) = Zarrsource()
RA._sourcetrait(::ZD.ZarrDataset) = Zarrsource()
