RA.sourcetrait(::NCD.Dataset) = NCDsource()
RA.sourcetrait(::NCD.Variable) = NCDsource()
RA.sourceconstructor(::NCDsource) = NCD.Dataset
RA.checkfilename(::NCDsource, filename) =
    isfile(filename) || RA._isurl(filename) || RA._filenotfound_error(filename)

# NCDatasets owns a file handle that must be released explicitly.
RA._close_dataset(ds::NCD.Dataset) = NCD.close(ds)
Base.close(os::RA.OpenStack{NCDsource}) = NCD.close(RA.dataset(os))

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