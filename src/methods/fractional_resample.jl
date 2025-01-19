using DimensionalData

fractional_resample(x, res; kw...) = fractional_resample(x; res, kw...)
fractional_resample(xs::RasterStackOrArray...; kw...) = fractional_resample(xs; kw...)
function fractional_resample(ser::AbstractRasterSeries, args...; kw...)
    map(x -> resample(x, args...; kw...), ser)
end
function fractional_resample(A::AbstractRaster; 
    categories::NamedTuple, to=nothing, res=nothing, size=nothing
)
    ds = _extent2dims(to; size, res, crs(A))
    shared_dims = commondims(A, ds)
    isintervals(shared_dims) || throw(ArgumentError("Fractional resampling only supported for `Intervals`."))
    outputdims = setdims(dims(A), ds)
    outputdata = Array{typeof(categories)}(undef, outputdims)
    output = rebuild(A; data=outputdata, dims=outputdims, refdims=())
    map(DimIndices(otherdims(A, ds))) do D
        slice = view(A, D...)
        fractions = read_fractions(slice, categories, weights)
    end
end

function read_fractions(slice, categories, weights, indices::CartesianIndices)
    vals = view(slice, indices)
    map(categories) do c
        sum(zip(vals, weights)) do (v, w)
            _compare(c, v) * w
        end / sum(weights)
    end
end

_compare(c::Interval, v) = v in c
_compare(c, v) = v == c
