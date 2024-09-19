function RA.OpenStack(fs::RA.FileStack{Zarrsource,K}) where K
    RA.OpenStack{Zarrsource,K}(ZD.ZarrDataset(RA.filename(fs)), fs.mods)
end

# In ZarrDatasets, the file is open for reading the values and closed afterwards. 
Base.close(os::RA.OpenStack{Zarrsource}) = nothing

function RA._open(f, ::Zarrsource, filename::AbstractString; write=false, kw...)
    ds = ZarrDatasets.ZarrDataset(filename)
    RA._open(f, Zarrsource(), ds; kw...)
end

RA._sourcetrait(::ZD.ZarrVariable) = Zarrsource()
RA._sourcetrait(::ZD.ZarrDataset) = Zarrsource()

RA.missingval(var::ZD.ZarrVariable, args...) = RA.missingval(RA.Metadata{Zarrsource}(CDM.attribs(var)))
RA.missingval(var::ZD.ZarrVariable, md::RA.Metadata{<:Zarrsource}) = RA.missingval(md)

@inline function RA.get_scale(metadata::RA.Metadata{<: Zarrsource}, scaled::Bool)
    scale = scaled ? get(metadata, "scale_factor", nothing) : nothing
    offset = scaled ? get(metadata, "add_offset", nothing) : nothing
    return scale, offset
end
# TODO: handle multiple missing values
function RA.missingval(md::RA.Metadata{<:Zarrsource})
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
