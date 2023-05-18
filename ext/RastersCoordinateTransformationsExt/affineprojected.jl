import Rasters: AffineProjected

function AffineProjected(f;
    data=AutoIndex(), metadata=NoMetadata(), crs=nothing, mappedcrs=nothing, paired_lookup, dim=AutoDim()
)
    AffineProjected(f, data, metadata, crs, mappedcrs, paired_lookup, dim)
end

DD.dim(lookup::AffineProjected) = lookup.dim
RA.crs(lookup::AffineProjected) = lookup.crs
RA.mappedcrs(lookup::AffineProjected) = lookup.mappedcrs
paired_lookup(lookup::AffineProjected) = lookup.paired_lookup

DD.metadata(lookup::AffineProjected) = lookup.metadata
function DD.rebuild(l::AffineProjected;
    affinemap=l.affinemap, data=l.data, metadata=metadata(l),
    crs=crs(l), mappedcrs=mappedcrs(l), paired_lookup=paired_lookup(l), dim=dim(l), args...
)
    AffineProjected(affinemap, data, metadata, crs, mappedcrs, paired_lookup, dim)
end
function Dimensions.format(l::AffineProjected, D::Type, index, axis::AbstractRange)
    return rebuild(l; data=axis, dim=basetypeof(D)())
end
DD.bounds(lookup::AffineProjected) = _bounds(dim(lookup), lookup)

function _bounds(::XDim, lookup::AffineProjected)
    am = lookup.affinemap
    extrema = _affine_extrema(am, lookup, lookup.paired_lookup)
    xbounds = min(map(first, extrema)...), max(map(first, extrema)...)
    return xbounds
end
function _bounds(::YDim, lookup::AffineProjected)
    am = lookup.affinemap
    extrema = _affine_extrema(am, lookup.paired_lookup, lookup)
    ybounds = min(map(last, extrema)...), max(map(last, extrema)...)
    return ybounds
end

LA.transformfunc(lookup::AffineProjected) = CoordinateTransformations.inv(lookup.affinemap)

function Dimensions.sliceunalligneddims(
    _, I::NTuple{2,<:Union{Colon,AbstractArray}},
    ud1::Dimension{<:AffineProjected}, ud2::Dimension{<:AffineProjected}
)
    # swap colons for the dimension index, which is the same as the array axis
    I = map((ud1, ud2), I) do d, i
        i isa Colon ? parent(lookup(d)) : i
    end

    af = lookup(ud1).affinemap
    # New resolution for step size changes
    M = af.linear * [step(I[1]) 0; 0 step(I[2])]
    # New of origin for slice
    v = af(collect(first.(I) .- 1))
    # Create a new affine map
    affinemap = CoordinateTransformations.AffineMap(M, v)
    # Build new lookups with the affine map. Probably should define `set` to do this.
    newlookups = map((ud1, ud2), I) do d, i
        rebuild(lookup(d); data=Base.OneTo(length(i)), affinemap)
    end
    newdims = map((ud1, ud2), newlookups, map(parent, reverse(newlookups))) do d, lookup, paired_lookup
        rebuild(d, rebuild(lookup; paired_lookup))
    end
    newrefdims = ()
    return newdims, newrefdims
end

function Base.reverse(lookup::AffineProjected)
    sp = reverse(span(lookup))
    rebuild(lookup; order=o, span=sp)
end


function Base.show(io::IO, mime::MIME"text/plain", lookup::AffineProjected)
    LA.show_compact(io, mime, lookup)
    print(io, " ")
    printstyled(io, lookup.affinemap; color=:cyan)
    if !(crs(lookup) isa Nothing)
        print(io, " crs: ", nameof(typeof(crs(lookup))))
    end
    if !(mappedcrs(lookup) isa Nothing)
        print(io, " mappedcrs: ", nameof(typeof(mappedcrs(lookup))))
    end
end

function _dims2geotransform(x::XDim{<:AffineProjected}, y::YDim)
    _affine2geotransform(parent(x).affinemap)
end

function Rasters.maybe_resample(l::AffineProjected, A)
    res = Y(abs(l.affinemap.linear[1, 1])), X(abs(l.affinemap.linear[2, 2]))
    A = resample(A, res)
end
