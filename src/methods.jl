"""
    replace_missing(a::AbstractGeoArray, newmissing)

Replace missing values in the array with a new missing value, also
updating the missingval field.
"""
replace_missing(a::AbstractGeoArray, newmissing=missing) = begin
    newdata = if ismissing(missingval(a))
        collect(Missings.replace(data(a), newmissing))
    else
        replace(data(a), missingval(a) => newmissing)
    end
    rebuild(a; data=newdata, missingval=newmissing)
end

# Helper methods ##############################################################
boolmask(A::AbstractArray) = boolmask(A, missing)
boolmask(A::AbstractGeoArray) =
    rebuild(A; data=boolmask(A, missingval(A)), missingval=false, name="Boolean mask")
boolmask(A::AbstractGeoArray, missingval::Missing) =
    (x -> !ismissing(x)).(data(A))
boolmask(A::AbstractGeoArray, missingval) =
    (x -> !isapprox(x, missingval)).(data(A))

missingmask(A::AbstractArray) = missingmask(A, missing)
missingmask(A::AbstractGeoArray) =
    rebuild(A; data=missingmask(A, missingval(A)), missingval=missing, name="Missing mask")
missingmask(A::AbstractGeoArray, missingval::Missing) =
    (a -> ismissing(a) ? missing : true).(data(A))
missingmask(A::AbstractGeoArray, missingval) =
    (a -> isapprox(a, missingval) ? missing : true).(data(A))
