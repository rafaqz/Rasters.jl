# Also see the DimentionalData.jl interface

"""
    missingval(x)

Returns the value representing missing data in the dataset
"""
function missingval end
missingval(x) = missing

"""
Get the crs projection of a dim or for the `Lat`/`Lon` dims of an array.
"""
function crs end

"""
    usercrs(x)

Get the user facing crs projection of a dim or for the `Lat`/`Lon` dims of an array.

This is used to convert `Selector` values form the user defined projection 
to the underlying projection, and to show plot axes in the user projection.
"""
function usercrs end

"""
    dimcrs(x)

Get the index crs projection of a dim or for the `Lat`/`Lon` dims of an array.

Where the dimension mode is `Converted`. This is often used in netcdf where
the underlying projection of the data is not what is contained in the vector index.
"""
