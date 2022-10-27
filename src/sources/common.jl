# Common utils methods for building Netcdf and GRIB Rasters.
const NCD = NCDatasets
const GDS = GRIBDatasets

const NCD_OR_GDS = Union{NCD.Dataset, GDS.Dataset}

_unuseddimerror(dimname) = error("Netcdf contains unused dimension $dimname")

function _dsdim(ds::NCD.Dataset, dimname::Key, crs=nothing, mappedcrs=nothing)
    if haskey(ds, dimname)
        D = _dsdimtype(dimname)
        lookup = _dslookup(ds, dimname, D, crs, mappedcrs)
        return D(lookup)
    else
        # The var doesn't exist. Maybe its `complex` or some other marker,
        # so make it a custom `Dim` with `NoLookup`
        len = _dsfinddimlen(ds, dimname)
        len === nothing && _unuseddimerror()
        return Dim{Symbol(dimname)}(NoLookup(Base.OneTo(len)))
    end
end

function _dsdim(ds::GDS.Dataset, dimname::Key, crs=nothing, mappedcrs=nothing)
    D = _gdsdimtype(dimname)
    lookup = _dslookup(ds, dimname, D, crs, mappedcrs)
    return D(lookup)
end

function _dsfinddimlen(ds, dimname)
    for key in keys(ds)
        var = NCD.variable(ds, key)
        dimnames = NCD.dimnames(var)
        if dimname in dimnames
            return size(var)[findfirst(==(dimname), dimnames)]
        end
    end
    return nothing
end

function _dsspan(index, order)
    # Handle a length 1 index
    length(index) == 1 && return Regular(zero(eltype(index)))
    step = index[2] - index[1]
    for i in 2:length(index)-1
        # If any step sizes don't match, its Irregular
        if !(index[i+1] - index[i] ≈ step)
            bounds = if length(index) > 1
                beginhalfcell = abs((index[2] - index[1]) * 0.5)
                endhalfcell = abs((index[end] - index[end-1]) * 0.5)
                if LA.isrev(order)
                    index[end] - endhalfcell, index[1] + beginhalfcell
                else
                    index[1] - beginhalfcell, index[end] + endhalfcell
                end
            else
                index[1], index[1]
            end
            return Irregular(bounds)
        end
    end
    # Otherwise regular
    return Regular(step)
end

# _dslookup
# Generate a `LookupArray` from a netcdf dim.
function _dslookup(ds::NCD_OR_GDS, dimname, D, crs, mappedcrs)
    dvar = ds[dimname]
    index = dvar[:]
    filetype = ds isa NCD.Dataset ? NCDfile : GRIBfile
    metadata = Metadata{filetype}(LA.metadatadict(dvar.attrib))
    return _dslookup(ds, dimname, D, index, metadata, crs, mappedcrs)
end
# For unknown types we just make a Categorical lookup
function _dslookup(ds::NCD_OR_GDS, dimname, D, index::AbstractArray, metadata, crs, mappedcrs)
    Categorical(index; metadata=metadata)
end
# For Number and AbstractTime we generate order/span/sampling
function _dslookup(
    ds::NCD.Dataset, dimname, D, index::AbstractArray{<:Union{Number,Dates.AbstractTime}},
    metadata, crs, mappedcrs
)
    # Assume the locus is at the center of the cell if boundaries aren't provided.
    # http://cfconventions.org/cf-conventions/cf-conventions.html#cell-boundaries
    order = _dsorder(index)
    var = NCD.variable(ds, dimname)
    if haskey(var.attrib, "bounds")
        boundskey = var.attrib["bounds"]
        boundsmatrix = Array(ds[boundskey])
        span, sampling = Explicit(boundsmatrix), Intervals(Center())
        return _dslookup(D, index, order, span, sampling, metadata, crs, mappedcrs)
    elseif eltype(index) <: Dates.AbstractTime
        span, sampling = _ncdperiod(index, metadata)
        return _dslookup(D, index, order, span, sampling, metadata, crs, mappedcrs)
    else
        span, sampling = _dsspan(index, order), Points()
        return _dslookup(D, index, order, span, sampling, metadata, crs, mappedcrs)
    end
end

function _dslookup(
    ds::GDS.Dataset, dimname, D, index::AbstractArray{<:Union{Number,Dates.AbstractTime}},
    metadata, crs, mappedcrs
)
    # Assume the locus is at the center of the cell if boundaries aren't provided.
    # http://cfconventions.org/cf-conventions/cf-conventions.html#cell-boundaries
    order = _dsorder(index)
    # var = NCD.variable(ds, dimname)
    if dimname in ["time", "valid_time"]
        # We consider the epoch 1970-01-01T00:00:00, as it appears to be in gribs files
        # dates = Dates.unix2datetime.(index)
        dates = Second.(index) .+ GDS.DEFAULT_EPOCH
        steps = unique(dates[2:end] .- dates[1:end-1])
        if length(steps) == 1
            span, sampling = Regular(steps[1]), Points()
        else
            span, sampling = Irregular((extrema(dates)...)), Points()
        end
        index = dates
        return _dslookup(D, index, order, span, sampling, metadata, crs, mappedcrs)
    else
        span, sampling = _dsspan(index, order), Points()
        return _dslookup(D, index, order, span, sampling, metadata, crs, mappedcrs)
    end
end

# For X and Y use a Mapped <: AbstractSampled lookup
function _dslookup(
    D::Type{<:Union{<:XDim,<:YDim}}, index, order::Order, span, sampling, metadata, crs, mappedcrs
)
    # If the index is regularly spaced and there is no crs
    # then there is probably just one crs - the mappedcrs
    crs = if crs isa Nothing && span isa Regular
        mappedcrs
    else
        crs
    end
    dim = DD.basetypeof(D)()
    return Mapped(index, order, span, sampling, metadata, crs, mappedcrs, dim)
end

# Band dims have a Categorical lookup, with order
function _dslookup(D::Type{<:Band}, index, order::Order, span, sampling, metadata, crs, mappedcrs)
    Categorical(index, order, metadata)
end
# Otherwise use a regular Sampled lookup
function _dslookup(D::Type, index, order::Order, span, sampling, metadata, crs, mappedcrs)
    Sampled(index, order, span, sampling, metadata)
end

function _dsorder(index)
    index[end] > index[1] ? ForwardOrdered() : ReverseOrdered()
end