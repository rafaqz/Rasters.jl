import Rasters: DD, Extents, GI
import Rasters: OBJ_ARGUMENT, GEOMETRYCOLUMN_KEYWORD, _get_geometries, _run, RasterStackOrArray
"""
    mapzonal(reducer, operator, data::RasterStack; of, kw...)

Think of this like `zonal`, but a named tuple gets sent to the `operator`, 
and `reducer` gets the result of `operator.(data[mask])`.

# Arguments

- `reducer`: any function that can consume an iterable as its only argument, such as `sum` or `Statistics.mean`.  Can also be `identity` or any other user-defined function!
- `operator`: any function that can consume a named tuple as its only argument, such as `mean` or `std`.  Applied on the set of values from RasterStacks, what you would get when you index into a RasterStack like `rasterstack[1, 1]`.
- `data`: a `RasterStack`
- `of`: a `DimTuple`, `Extent`, $OBJ_ARGUMENT

# Keywords
$GEOMETRYCOLUMN_KEYWORD
These can be used when `of` is or contains (a) GeoInterface.jl compatible object(s):

- `shape`: Force `data` to be treated as `:polygon`, `:line` or `:point`, where possible.
- `boundary`: for polygons, include pixels where the `:center` is inside the polygon,
    where the line `:touches` the pixel, or that are completely `:inside` inside the polygon.
    The default is `:center`.
- `progress`: show a progress bar, `true` by default, `false` to hide..
- `skipmissing`: wether to apply `f` to the result of `skipmissing(A)` or not. If `true`
    `f` will be passed an iterator over the values, which loses all spatial information.
    if `false` `f` will be passed a masked `Raster` or `RasterStack`, and will be responsible
    for handling missing values itself. The default value is `true`.

# Example

"""
function mapzonal(reducer, operator, data::RasterStackOrArray; of, kw...)
    _mapzonal(reducer, operator, data, of; kw...)
end

function _mapzonal(reducer, operator, x, ext::Extents.Extent; skipmissing = true)
    cropped = crop(x; to = ext, touches = true)
    prod(size(cropped)) > 0 || return missing
    # We can't use skipmissing here, since it doesn't work on rasterstacks
    if skipmissing
        return reducer( # reduce the result of the following operations - many reducers don't support iterators.
            map( # apply operator to each value
                operator, 
                Iterators.filter( # skip missing values
                    Base.Fix1(any, !ismissing), 
                    Iterators.map( # get the values as named tuples from the rasterstack
                        Base.Fix1(getindex, cropped), 
                        DD.DimIndices(cropped)
                    )
                )
            )
        )
    else
        return reducer( # apply the reducer function
            Iterators.map(DD.DimIndices(cropped)) do I
                operator(cropped[I]) # get the value as namedtuples from the rasterstack and apply the operator
            end
        )
    end
end

function _mapzonal(reducer, operator, x, of; kw...)
    # Otherwise of is a geom, table or vector
    _mapzonal(reducer, operator, x, GI.trait(of), of; kw...)
end

function _mapzonal(reducer, operator, x, ::GI.AbstractFeatureCollectionTrait, fc; kw...)
    _mapzonal(reducer, operator, x, nothing, fc; kw...) # treat this as a table of geometries
end

# This handles tables, featurecollections and vectors of geometries.
function _mapzonal(reducer, operator, x, ::Nothing, data; progress = true, threaded = true, geometrycolumn = nothing, kw...)
    geoms = _get_geometries(data, geometrycolumn)
    n = length(geoms)
    n == 0 && return []
    zs, start_index = _alloc_mapzonal(reducer, operator, x, geoms, n; kw...)
    _run(start_index:n, threaded, progress, "Applying $reducer and $operator to each geometry...") do i 
        zs[i] = _mapzonal(reducer, operator, x, geoms[i]; kw...)
    end
    return zs
end

function _alloc_mapzonal(reducer, operator, x, geoms, n; kw...)
    # Find first non-missing entry and count number of missing entries
    n_missing::Int = 0
    z1 = _mapzonal(reducer, operator, x, first(geoms); kw...)
    for geom in geoms
        z1 = _mapzonal(reducer, operator, x, geom; kw...)
        if !ismissing(z1)
            break
        end
        n_missing += 1
    end
    zs = Vector{Union{Missing,typeof(z1)}}(undef, n)
    zs[1:n_missing] .= missing
    # Exit early when all elements are missing
    if n_missing == n
        return zs, n_missing + 1
    end
    zs[n_missing + 1] = z1
    return zs, n_missing + 1
end

# This handles single features (just decomposes to geometry)
function _mapzonal(reducer, operator, data, ::GI.AbstractFeatureTrait, feature; kw...)
    _mapzonal(reducer, operator, data, GI.geometry(feature); kw...)
end
    

# Now, we get into the meat of handling actual geometry

function _mapzonal(reducer, operator, st::RasterStackOrArray, ::GI.AbstractGeometryTrait, geom; 
    skipmissing=true, kw...
)
    cropped = crop(st; to=geom, touches=true)
    prod(size(cropped)) > 0 || return missing # mapzonal should always return ONE value...
    # Construct a "boolean mask" of the same size as `cropped`, 
    # with `true` wherever the pixel is inside `geom`.
    # This is controlled by the `rasterize` keywords that `boolmask` accepts.
    mask_raster = boolmask(geom; to=cropped, kw...)
    indices = view(DD.DimIndices(mask_raster), mask_raster)
    # We can now use this pre-made boolean mask to index into `cropped`.
    if skipmissing # TODO: decide on whether to use map or Iterators.map.  Iterators version is faster and allocates less, but is less generally applicable.
        # Maybe users can do `sum \circ collect` if they want to get a vector??
        return reducer(map(operator, Iterators.filter(Base.Fix1(any, !ismissing), Iterators.map(Base.Fix1(getindex, cropped), indices))))
    else
        return reducer(map(view(cropped, mask)) do val
            operator(val)
        end)
    end
end