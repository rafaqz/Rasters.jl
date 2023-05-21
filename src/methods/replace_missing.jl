"""
    replace_missing(a::AbstractRaster, newmissingval)
    replace_missing(a::AbstractRasterStack, newmissingval)

Replace missing values in the array or stack with a new missing value,
also updating the `missingval` field/s.

# Keywords

- `filename`: a filename to write to directly, useful for large files.
- `suffix`: a string or value to append to the filename.
    A tuple of `suffix` will be applied to stack layers. `keys(st)` are the default.

# Example

```jldoctest
using Rasters, RasterDataSources
A = Raster(WorldClim{Climate}, :prec; month=1) |> replace_missing
missingval(A)
# output
missing
```

"""
replace_missing(x; missingval=missing, kw...) = replace_missing(x, missingval; kw...)
function replace_missing(A::AbstractRaster{T}, missingval::MV;
    filename=nothing, suffix=nothing
) where {T,MV}
    MT = if ismissing(missingval)
        promote_type(T, Missing)
    else
        promote_type(nonmissingtype(T), MV)
    end
    old_missingval = Rasters.missingval(A)
    missingval = convert(MT, missingval)
    repmissing(x) = isequal(x, old_missingval) || ismissing(x) ? missingval : x
    # Disk-backed arrays need to be lazy, memory-backed don't.
    # But in both cases we make sure we return an array with the missingval
    # in the eltype, even if there are no missing values in the array.
    if !isnothing(filename)
        A1 = create(filename, MT, dims(A); 
            parent=parent(A), suffix, missingval, name=name(A), metadata=metadata(A)
        )
        open(A1; write=true) do O
            O .= repmissing.(A)
        end
        return A1
    else
        # We need to force T of Union{T,Missing} for DiskArrays broadcasts
        if isdisk(A)
            data = repmissing.(parent(A))
            if missingval isa Missing
                data = MissingDiskArray(MT, data)
            end
        else
            data = similar(parent(A), MT)
            data .= repmissing.(parent(A))
        end
        return rebuild(A; data, missingval)
    end
end
function replace_missing(st::AbstractRasterStack, args...; suffix=keys(st), kw...)
    mapargs(st, suffix) do A, s
        replace_missing(A, args...; suffix=s, kw...)
    end
end
function replace_missing(s::AbstractRasterSeries, args...; kw...)
    map(x -> replace_missing(x, args...; kw...), s)
end

function replace_missing!(s::RasterSeriesOrStack, args...; kw...)
    map(x -> replace_missing!(x, args...; kw...), s)
end
function replace_missing!(A::AbstractRaster, missingval=missing)
    repmissing(x) = isequal(x, Rasters.missingval(A)) ? missingval : x
    A .= repmissing.(A)
    return rebuild(A; missingval)
end

