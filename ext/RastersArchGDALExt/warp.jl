function warp(A::AbstractRaster, flags::Dict; filename=nothing, kw...)
    odims = otherdims(A, (X, Y, Band))
    if length(odims) > 0
        isnothing(filename) || throw(ArgumentError("Cannot currently write dimensions other than X/Y/Band to disk using `filename` keyword. Make a Rasters.jl github issue if you need this."))
        # Handle dimensions other than X, Y, Band
        slices = slice(A, odims)
        warped = map(A -> _warp(A, flags; kw...), slices)
        return combine(warped, odims)
    else
        return _warp(A, flags; filename, kw...)
    end
end
function warp(st::AbstractRasterStack, flags::Dict; filename=nothing, suffix=keys(st), kw...)
    RA.mapargs((A, s) -> warp(A, flags; filename, suffix=s), st, suffix; kw...)
end

function _warp(A::AbstractRaster, flags::Dict; 
    filename=nothing, 
    suffix="", 
    missingval=nokw,
    maskingval=Rasters.missingval(A),
    name=Rasters.name(A),
    kw...
)
    A1 = _set_gdalwarp_sampling(A)
    filename = RA._maybe_add_suffix(filename, suffix)
    flagvect = reduce([flags...]; init=String[]) do acc, (key, val)
        append!(acc, String[_asflag(key), _stringvect(val)...])
    end
    # TODO: detect if `A` already holds a lazy GDAL FileArray. 
    # If it does, we can just open it and use it directly.
    tempfile = isnothing(filename) ? nothing : tempname() * ".tif"
    warp_kw = isnothing(filename) || filename == "/vsimem/tmp" ? () : (; dest=filename)
    # We really need a missingval for `warp`, as it may rotate and add missing value
    missingval = if RA.isnokw(missingval) 
        if RA.missingval(A) isa Union{Missing,Nothing} 
            RA._type_missingval(Missings.nonmissingtype(eltype(A)))
        else
            RA.missingval(A)
        end
    else
        missingval
    end
    out = AG.Dataset(A1; filename=tempfile, missingval, kw...) do dataset
        AG.gdalwarp([dataset], flagvect; warp_kw...) do warped
            # Read the raster lazily, dropping Band if there is none in `A`
            raster = Raster(warped; lazy=true, dropband=!hasdim(A, Band()), name, maskingval)
            # Either read the MEM dataset to an Array, or keep a filename base raster lazy
            return isnothing(filename) ? read(raster) : raster
        end
    end
    # And permute the dimensions back to what they were in A
    out1 = _maybe_restore_from_gdal(out, dims(A))
    out2 = _reset_gdalwarp_sampling(out1, A)
    return out2
end

_asflag(x) = string(x)[1] == '-' ? x : string("-", x)

_stringvect(x::AbstractVector) = Vector(string.(x))
_stringvect(x::Tuple) = [map(string, x)...]
_stringvect(x) = [string(x)]

function _set_gdalwarp_sampling(A)
    x = if sampling(A, X) isa Points
        DD.maybeshiftlocus(Start(), set(convertlookup(Projected, dims(A, X)), Intervals(Center())))
    else
        DD.maybeshiftlocus(Start(), convertlookup(Projected, dims(A, X)))
    end
    y = if sampling(A, Y) isa Points
        DD.maybeshiftlocus(Start(), set(convertlookup(Projected, dims(A, Y)), Intervals(Center())))
    else
        DD.maybeshiftlocus(Start(), convertlookup(Projected, dims(A, Y)))
    end
    return set(A, X => x, Y=> y)
end

function _reset_gdalwarp_sampling(A, template)
    x = if sampling(template, X) isa Points
        set(DD.maybeshiftlocus(Center(), lookup(A, X)), Points())
    else
        DD.maybeshiftlocus(locus(template, X), lookup(A, X))
    end
    y = if sampling(template, Y) isa Points
        set(DD.maybeshiftlocus(Center(), lookup(A, Y)), Points())
    else
        DD.maybeshiftlocus(locus(template, Y), lookup(A, Y))
    end
    return set(A, X => x, Y => y)
end
