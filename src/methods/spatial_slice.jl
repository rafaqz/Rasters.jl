# Spatial-slice helpers shared by `zonal` and `extract`.
#
# These functions handle the "apply `f` to each spatial slice of a higher-rank
# raster" pattern: when the raster has dims beyond X/Y, we want to either map
# `f` over each X/Y slice (and return a vector data cube), or concatenate
# per-geometry sliced results into a single backing array.

# Apply `f` to each spatial slice of `x`. The non-spatial dims are kept; the
# spatial dims are reduced inside `f`. This differs from `mapslices` in that
# the sliced dims are dropped automatically, which is what we want.
# `AbstractDimArray` is required only because we use `DD.otherdims`.
function _mapspatialslices(f, x::AbstractDimArray;
    spatialdims=(Val{DD.XDim}(), Val{DD.YDim}()),
    missingval=missingval(x),
)
    dimswewant = DD.otherdims(x, spatialdims)
    if isempty(dimswewant)
        return f(x)
    end
    slicedims = rebuild.(dims(x, dimswewant), axes.((x,), dimswewant))
    if any(isempty, DD.dims(x, spatialdims))
        # If any of the spatial dims are empty, we can just return a constant
        # missing array - this way we don't construct the dimslices at all.
        missing_array = FillArrays.Fill{Union{typeof(missingval),eltype(x)},length(dimswewant)}(
            missingval, length.(dimswewant)
        )
        return rebuild(x; data=missing_array, dims=dimswewant, refdims=())
    end
    iterator = (rebuild(x; data=d, dims=dims(d)) for d in DD.DimSlices(x; dims=slicedims, drop=true))
    return rebuild(x; data=f.(iterator), dims=dimswewant, refdims=())
end
# SkipMissingVal and SkipMissing both store the initial value in `x`.
function _mapspatialslices(f, s::Union{SkipMissingVal,Base.SkipMissing};
    spatialdims=(Val{DD.XDim}(), Val{DD.YDim}()),
    missingval=missingval(s.x),
)
    return _mapspatialslices(f ∘ skipmissing, s.x; spatialdims, missingval)
end

# Wrap `f` so it applies to each spatial slice. Returns `f` unchanged when
# `spatialslices` is false.
_maybe_spatialsliceify(f, spatialslices, missingval=missing) =
    istrue(spatialslices) ? _SpatialSliceify(f, (Val{DD.XDim}(), Val{DD.YDim}()), missingval) : f
_maybe_spatialsliceify(f, spatialslices::DD.AllDims, missingval=missing) =
    _SpatialSliceify(f, spatialslices, missingval)

"""
    _SpatialSliceify(f, dims)

A callable struct that applies `mapslices(f, x; dims = spatialdims)` to the
input array `x`, and removes empty dimensions.

```jldoctest
data = ones(10, 10, 10, 10);
f = _SpatialSliceify(sum, (1, 2))
size(f(data))

# output
(10, 10)
```
"""
struct _SpatialSliceify{F,D,M}
    f::F
    dims::D
    missingval::M
end

(r::_SpatialSliceify)(x) = _mapspatialslices(r.f, x; spatialdims=r.dims, missingval=r.missingval)

# Concatenate an array of arrays along a new last dimension. Returns a
# `ConcatDiskArray` so we don't allocate; users can rechunk later.
function __do_cat_with_last_dim(input_arrays)
    # Assumes the input is a vector of arrays.
    As = Missings.disallowmissing(collect(input_arrays))
    dims = ndims(first(As)) + 1
    sz = ntuple(dims) do i
        i == dims ? length(As) : 1
    end
    cdas = reshape(As, sz)
    return DiskArrays.ConcatDiskArray(cdas)
end

function __do_cat_with_last_dim_multidim_version(As)
    # Does NOT assume the input is a vector of arrays.
    new_n_dims = ndims(As) + ndims(first(As))
    sz = ntuple(new_n_dims) do i
        i <= ndims(first(As)) ? 1 : size(As, i - 1)
    end
    cdas = reshape(As, sz)
    return DiskArrays.ConcatDiskArray(cdas)
end

# Cat children along a new dim and rebuild as a dim array off `parent`.
function _cat_and_rebuild_parent(parent, children, newdim)
    backing_array = __do_cat_with_last_dim(children)
    children_dims = dims(first(children))
    final_dims = DD.format((children_dims..., newdim), backing_array)
    return rebuild(parent; data=backing_array, dims=final_dims)
end

precompile(__do_cat_with_last_dim, (Vector{Raster{<:Any,1}},))
precompile(__do_cat_with_last_dim, (Vector{Raster{<:Any,2}},))
precompile(__do_cat_with_last_dim, (Vector{Raster{<:Any,3}},))
precompile(__do_cat_with_last_dim, (Vector{Raster{<:Any,4}},))
