"""
    read(A::AbstractGeoArray)

`read` will move any `AbstractGeoArray` completely to memory, as a `GeoArray`.
```
"""
Base.read(A::AbstractGeoArray) = GeoArray(A)

"""
    read(A::AbstractGeoStack)

`read` will move any `AbstractGeoStack` completely to memory, as a `GeoStack` of `GeoArray`.
"""
Base.read(st::AbstractGeoStack) = GeoStack(st)

"""
    read(A::AbstractGeoSeries)

`read` will move any `AbstractGeoSeries` completely to memory, as a `GeoSeries` 
of `GeoStack` or `GeoArray`.
"""
Base.read(ser::AbstractGeoSeries) = rebuild(ser; data=[read(ser[i]) for i in eachindex(ser)])

