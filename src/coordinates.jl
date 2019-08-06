coords(a, dims...) = coords(dims(a), sortdims(a, dims))
coords(coorddims::Tuple, indexdims) = 
    (coords(coorddims[1], indexdims[1])..., coords(tail(coorddims), tail(indexdims))...)
coords(coorddims::Tuple{}, indexdims::Tuple{}) = ()

coords(coorddim::AbstractDimension, indexdim::Nothing) = ()
coords(coorddim::AbstractDimension, indexdim::AbstractDimension) = 
    (basetype(coorddim)(collect(val(coorddim))[val(indexdim)]),)


# Reference rom GeoArrays
# Generate upper left coordinates for specic index
# coords(aff::AbstractAffineDims, point) = val(aff)(point .- 1)

# Generate center coordinates for specific index
# centercoords(aff::AbstractAffineDims, point) = val(aff)(point .- 0.5)

# Convert coordinates to indices
# indices(aff::AbsractAffineDims, point) = map(x -> round(Int, x), inv(val(aff))(point)) .+ 1
