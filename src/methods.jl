"""
    replace_missing(a::AbstractGeoArray, newmissingval)
    replace_missing(a::AbstractGeoStack, newmissingval)

Replace missing values in the array or stack with a new missing value, 
also updating the `missingval` field/s.
"""
replace_missing(a::DiskGeoArray, args...) = 
    replace_missing(GeoArray(a), args...)
replace_missing(a::MemGeoArray, newmissingval=missing) = begin
    newdata = if ismissing(missingval(a))
        collect(Missings.replace(parent(a), newmissingval))
    else
        replace(parent(a), missingval(a) => newmissingval)
    end
    rebuild(a; data=newdata, missingval=newmissingval)
end
replace_missing(stack::AbstractGeoStack, newmissingval=missing) = 
    rebuild(stack, map(a -> replace_missing(a, newmissingval, values(stack))))

"""
    boolmask(A::AbstractArray, [missingval])

Create a mask array of `Bool` values, from any AbstractArray. For `AbstractGeoArray` 
the default `missingval` is `missingval(A)`, for all other `AbstractArray`s 
it is `missing`.

The array returned from calling `boolmask` on a `AbstractGeoArray` is a 
[`GeoArray`](@ref) with the same size and fields as the oridingl array
"""
function boolmask end
boolmask(A::AbstractGeoArray) =
    rebuild(A; data=boolmask(A, missingval(A)), missingval=false, name=:Bool_mask)
boolmask(A::AbstractArray, missingval::Missing=missing) = (a -> !ismissing(a)).(parent(A))
boolmask(A::AbstractArray, missingval) =
    if isnan(missingval)
        (a -> !isnan(a)).(parent(A))
    else
        (a -> !isapprox(a, missingval)).(parent(A))
    end

"""
    missingmask(A::AbstractArray, [missingval])

Create a mask array of `missing` or `true` values, from any AbstractArray. 
For `AbstractGeoArray` the default `missingval` is `missingval(A)`, 
for all other `AbstractArray`s it is `missing`.

The array returned from calling `boolmask` on a `AbstractGeoArray` is a 
[`GeoArray`](@ref) with the same size and fields as the oridingl array
"""
function missingmask end
missingmask(A::AbstractGeoArray) =
    rebuild(A; data=missingmask(A, missingval(A)), missingval=missing, name=:Missing_mask)
missingmask(A::AbstractArray, missingval::Missing=missing) =
    (a -> ismissing(a) ? missing : true).(parent(A))
missingmask(A::AbstractArray, missingval) =
    if isnan(missingval)
        (a -> isnan(a) ? missing : true).(parent(A))
    else
        (a -> isapprox(a, missingval) ? missing : true).(parent(A))
    end
