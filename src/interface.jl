# Also see the DimentionalData.jl interface

"""
    missingval(x)

Returns the value representing missing data in the dataset
"""
function missingval end
missingval(x) = missing

"""
    safeapply(f::Function, ::AbstractGeoStack, source)

Wrapper method to apply a function to data object provided for by a data source.

This facilitates wrapping the custom file open/close requirements of specific 
source libraries to safely deal with disk or api sourced datasets.
"""
function safeapply end

"""
Get the crs projection of a dim or for the `Lat`/`Lon` dims of an array.
"""
function crs end

"""
Get the user facing crs projection of a dim or for the `Lat`/`Lon` dims of an array.

This is used to convert `Selector` values form the user defined projection 
to the underlying projection, and to show plot axes in the user projection.
"""
function usercrs end

function childtype end
