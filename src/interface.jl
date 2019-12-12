# Also see the DimentionalData.jl interface

"""
    missingval(x)
Returns the value representing missing data in the dataset
"""
function missingval end
missingval(x) = missing

"""
Replace missing values
"""
function replace_missing end

"""
    safeapply(f::Function, ::AbstractGeoStack, source)

Wrapper method to apply a function to data object provided for by a data source.

This facilitates wrapping the custom file open/close requirements of specific 
source libraries to safely deal with disk or api sourced datasets.
"""
function safeapply end
