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
