RA.sourceconstructor(::Type{NCDsource}) = NCD.Dataset

function RA.checkmode(::NCDsource, filename, append::Bool, force::Bool)
    if append
        isfile(filename) ? "a" : "c"
    else
        RA.check_can_write(filename, force)
        "c"
    end
end

Base.close(os::RA.OpenStack{NCDsource}) = NCD.close(RA.dataset(os))

function RA._open(f, ::NCDsource, filename::AbstractString; write=false, kw...)
    isfile(filename) || RA._isurl(filename) || RA._filenotfound_error(filename)
    mode = write ? "a" : "r"
    NCD.Dataset(filename, mode) do ds
        RA._open(f, NCDsource(), ds; kw...)
    end
end

RA._sourcetrait(::NCD.Dataset) = NCDsource()
RA._sourcetrait(::NCD.Variable) = NCDsource()

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