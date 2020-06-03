"""
    replace_missing(a::AbstractGeoArray, newmissing)

Replace missing values in the array with a new missing value, also
updating the `missingval` field.
"""
replace_missing(a::DiskGeoArray, args...) = 
    replace_missing(GeoArray(a), args...)
replace_missing(a::MemGeoArray, newmissing=missing) = begin
    newdata = if ismissing(missingval(a))
        collect(Missings.replace(data(a), newmissing))
    else
        replace(data(a), missingval(a) => newmissing)
    end
    rebuild(a; data=newdata, missingval=newmissing)
end

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
    rebuild(A; data=boolmask(A, missingval(A)), missingval=false, name="Boolean mask")
boolmask(A::AbstractArray, missingval::Missing=missing) =
    (x -> !ismissing(x)).(data(A))
boolmask(A::AbstractArray, missingval) =
    (x -> !isapprox(x, missingval)).(data(A))

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
    rebuild(A; data=missingmask(A, missingval(A)), missingval=missing, name="Missing mask")
missingmask(A::AbstractArray, missingval::Missing=missing) =
    (a -> ismissing(a) ? missing : true).(data(A))
missingmask(A::AbstractArray, missingval) =
    (a -> isapprox(a, missingval) ? missing : true).(data(A))
