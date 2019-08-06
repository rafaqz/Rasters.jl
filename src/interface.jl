# At this stage only small extsions to the DimentionalData interface

"""
Returns the value representing missing data in the dataset
"""
function missingval end
missingval(x) = missing

"""
Returns metadata related to the dataset
"""
function metadata end
metatdata(x) = Dict()

"""
Replace missing values
"""
function replace_missing end

replace_missing(a::AbstractGeoArray, x) = begin
    newa = replace(a, missing=>x)
    rebuild(newa, parent(newa), dims(newa), refdims(newa), missingval(newa))
end
