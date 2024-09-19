resample(x, res; kw...) = resample(x; res, kw...)
resample(xs::RasterStackOrArray...; kw...) = resample(xs; kw...)
function resample(ser::AbstractRasterSeries, args...; kw...)
    map(x -> resample(x, args...; kw...), ser)
end
function resample(xs::Union{Tuple,NamedTuple}; to=first(xs), kw...)
    map(x -> resample(x; to, kw...), xs)
end
function resample(A::RasterStackOrArray; 
    to=nothing, 
    res=nothing, 
    crs=nothing, 
    size=nothing, 
    method=:near, 
    kw...
)
    (isnothing(size) || isnothing(res)) || _size_and_res_error()

    # Flags to send to `warp`, then to GDAL
    flags = Dict{Symbol,Any}()

    # Method
    flags[:r] = method

    # check if only to has been provided before overwriting arguments
    onlyto = !isnothing(to) && !isnothing(dims(to)) && !(to isa Extents.Extent) && isnothing(res) && isnothing(size) && isnothing(crs)

    # Extent
    if to isa Extents.Extent || isnothing(to) || isnothing(dims(to))
        to = isnothing(to) || to isa Extents.Extent ? to : GeoInterface.extent(to)
        if !isnothing(to)
            # Get the extent of geometries
            (xmin, xmax), (ymin, ymax) = to[(:X, :Y)]
            flags[:te] = [xmin, ymin, xmax, ymax]
        end
    else
        all(hasdim(to, (XDim, YDim))) || throw(ArgumentError("`to` must have both `XDim` and `YDim` dimensions to resample with GDAL"))

        # Set res from `to` if it was not already set
        if isnothing(res) && isnothing(size)
            todims = dims(to, (XDim, YDim))
            isregular(todims) || throw(ArgumentError("`to` has irregular dimensions. Provide regular dimensions, or explicitly provide `res` or `size`."))
            ysize, xsize = length.(todims)
            flags[:ts] = [ysize, xsize]
        end
        (xmin, xmax), (ymin, ymax) = bounds(to, (XDim, YDim))
        flags[:te] = [xmin, ymin, xmax, ymax]
    end

    # CRS
    crs = if isnothing(crs) 
        if to isa Extents.Extent
            nothing
        else
            # get crs from `to` or `x` if none was passed in
            isnothing(Rasters.crs(to)) ? Rasters.crs(A) : Rasters.crs(to)
        end
    else # issomething(crs)
        if crs isa String
            error("""
                Strings as CRS aren't yet supported.  
                Please pass a `GeoFormatTypes.jl` CRS format, like `ESRIWellKnownText`, `ProjString`, or similar.

                You can find out more about GeoFormatTypes.jl at https://juliageo.org/GeoFormatTypes.jl/stable/.
                """
            )
        end
        crs
    end
    if !isnothing(crs)
        wkt = convert(String, convert(WellKnownText, crs))
        flags[:t_srs] = wkt
        if isnothing(Rasters.crs(A))
            @warn "You have set a crs to resample to, but the object does not have crs so GDAL will assume it is already in the target crs. Use `newraster = setcrs(raster, somecrs)` to fix this."
        end
    end

    # Resolution
    if !isnothing(res)
        xres, yres = if res isa Real
            res, res
        elseif res isa Tuple{<:Dimension{<:Real},<:Dimension{<:Real}}
            map(val, dims(res, (YDim, XDim)))
        elseif res isa Tuple{<:Real,<:Real}
            reverse(res)
        else
            throw(ArgumentError("`res` must be a `Real`, or a 2 `Tuple` of `Real` or `Dimension`s wrapping `Real`. Got $res"))
        end

        flags[:tr] = [yres, xres]
    end

    # Size
    if !isnothing(size)
        xsize, ysize = if size isa Int
            size, size
        elseif size isa Tuple{<:Dimension{Int},<:Dimension{Int}}
            map(val, dims(size, (YDim, XDim)))
        elseif size isa Tuple{Int,Int}
            dimnum(A, XDim) > dimnum(A, YDim) ? size : reverse(size)
        else
            throw(ArgumentError("`size` must be a `Int`, or a 2 `Tuple` of `Int` or `Dimension`s wrapping `Int`. Got $size"))
        end
        flags[:ts] = [ysize, xsize]
    end

    # resample with `warp`
    resampled = warp(A, flags; kw...)

    # Return crs to the original type, from GDAL it will always be WellKnownText
    if !isnothing(crs)
        resampled = setcrs(resampled, crs)
    end

    # if only to is provided and it has dims, make sure dims are the exact same 
    if onlyto
        newdims = (commondims(to, XDim, YDim)..., otherdims(A, (XDim, YDim))...)
        resampled = rebuild(resampled; dims =newdims)
    end

    return resampled
end

_size_and_res_error() = throw(ArgumentError("Include only `size` or `res` keywords, not both"))
