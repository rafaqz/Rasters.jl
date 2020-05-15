using .NCDatasets

export NCDarray, NCDstack, NCDstackMetadata, NCDarrayMetadata, NCDdimMetadata

struct NCDstackMetadata{K,V} <: ArrayMetadata{K,V}
    val::Dict{K,V}
end

struct NCDarrayMetadata{K,V} <: ArrayMetadata{K,V}
    val::Dict{K,V}
end

struct NCDdimMetadata{K,V} <: DimMetadata{K,V}
    val::Dict{K,V}
end


# Utils ########################################################################

ncread(f, path::String) = NCDatasets.Dataset(f, path)

nondimkeys(dataset) = begin
    dimkeys = keys(dataset.dim)
    removekeys = if "bnds" in dimkeys
        dimkeys = setdiff(dimkeys, ("bnds",))
        boundskeys = map(k -> dataset[k].attrib["bounds"], dimkeys)
        union(dimkeys, boundskeys)
    else
        dimkeys
    end
    setdiff(keys(dataset), removekeys)
end

# Add a var array to a dataset before writing it.
ncwritevar!(dataset, A::AbstractGeoArray{T,N}) where {T,N} = begin
    A = reorderindex(A, Forward()) |>
        a -> reorderrelation(a, Forward())
    if ismissing(missingval(A))
        # TODO default _FillValue for Int?
        fillvalue = get(metadata(A), "_FillValue", NaN)
        A = replace_missing(A, convert(T, fillvalue))
    end
    # Define required dim vars
    for dim in dims(A)
        key = lowercase(name(dim))
        haskey(dataset.dim, key) && continue

        index = ncshiftindex(dim)

        if dim isa Lat || dim isa Lon
            dim = convertmode(Converted, dim)
        end
        md = metadata(dim)
        attribvec = [] #md isa Nothing ? [] : [val(md)...]
        defDim(dataset, key, length(index))
        println("writing key: ", key, " of type: ", eltype(index))
        defVar(dataset, key, index, (key,); attrib=attribvec)
    end
    # TODO actually convert the metadata type
    attrib = if metadata isa NCDarrayMetadata
        deepcopy(val(metadata(A)))
    else
        Dict()
    end
    # Remove stack metdata if it is attached
    pop!(attrib, "dataset", nothing)
    # Set missing value
    if !ismissing(missingval(A))
        attrib["_FillValue"] = convert(T, missingval(A))
    end
    key = name(A)
    println("writing key: ", key, " of type: ", T)
    dimnames = lowercase.(name.(dims(A)))
    attribvec = [attrib...]
    var = defVar(dataset, key, eltype(A), dimnames; attrib=attribvec)

    var[:] = data(A)

end

ncshiftindex(dim::Dimension) = ncshiftindex(mode(dim), dim)
ncshiftindex(mode::Sampled, dim::Dimension) = begin
    if span(mode) isa Regular
        if dim isa TimeDim
            if eltype(dim) isa Dates.AbstractDateTime
                val(dim)
            else
                shiftindexloci(dim, Start())
            end
        else
            shiftindexloci(dim, Center())
        end
    else
        val(dim)
    end
end
ncshiftindex(mode::IndexMode, dim::Dimension) = val(dim)

# CF standards don't enforce dimension names.
# But these are common, and should take care of most dims.
const dimmap = Dict("lat" => Lat,
                    "latitude" => Lat,
                    "lon" => Lon,
                    "long" => Lon,
                    "longitude" => Lon,
                    "time" => Ti,
                    "lev" => Vert,
                    "level" => Vert,
                    "vertical" => Vert,
                    "x" => X,
                    "y" => Y,
                    "z" => Z,
                   )

# Array ########################################################################
"""
    NCDarray(filename::AbstractString; name="", refdims=())

Create an array from a path to a netcdf file. The first non-dimension
layer of the file will be used as the array.

This is an incomplete implementation of the NetCDF standard. It will currently
only handle simple files in lattitude/longitude format. Real projections are
not yet handled.

If you need to use crs with NetCDF, make a fewture request in the issue queue.

## Arguments
- `filename`: `String` pointing to a netcdf file.

## Keyword arguments
- `name`: Name for the array. Will use array key if not supplied.
- `refdims`: Add dimension position array was sliced from. Mostly used programatically.
  loading the array. Can save on disk load time for large files.
"""
struct NCDarray{T,N,A,D<:Tuple,R<:Tuple,Na<:AbstractString,Me,Mi,S,K
               } <: DiskGeoArray{T,N,D,LazyArray{T,N}}
    filename::A
    dims::D
    refdims::R
    name::Na
    metadata::Me
    missingval::Mi
    size::S
    key::K
end
NCDarray(filename::AbstractString, key...; kwargs...) =
    ncread(dataset -> NCDarray(dataset, filename, key...; kwargs...), filename)
NCDarray(dataset::NCDatasets.Dataset, filename, key=nothing;
         refdims=(),
         dims=nothing,
         name=nothing,
         metadata=nothing,
         usercrs=nothing,
        ) = begin
    keys_ = nondimkeys(dataset)
    key = key isa Nothing || !(string(key) in keys_) ? first(keys_) : string(key)
    var = dataset[key]
    dims = dims isa Nothing ? GeoData.dims(dataset, key, usercrs) : dims
    name = name isa Nothing ? string(key) : name
    metadata_ = metadata isa Nothing ? GeoData.metadata(var, GeoData.metadata(dataset)) : metadata
    missingval = missing
    size_ = map(length, dims)
    T = eltype(var)
    N = length(dims)

    NCDarray{T,N,typeof.((filename,dims,refdims,name,metadata_,missingval,size_,key))...
       }(filename, dims, refdims, name, metadata_, missingval, size_, key)
end

key(A::NCDarray) = A.key

# AbstractGeoArray methods

withsourcedata(f, A::NCDarray) =
    ncread(dataset -> f(dataset[string(key(A))]), filename(A))

# Base methods

Base.size(A::NCDarray) = A.size

"""
    Base.write(filename::AbstractString, ::Type{NCDarray}, s::AbstractGeoArray)

Write an NCDarray to a netcdf file using NCDatasets.jl
"""
Base.write(filename::AbstractString, ::Type, A::AbstractGeoArray) = begin
    meta = metadata(A)
    if meta isa Nothing
        dataset = NCDatasets.Dataset(filename, "c")
    else
        # Remove the dataset metadata
        stackmd = pop!(deepcopy(val(meta)), "dataset", Dict())
        dataset = NCDatasets.Dataset(filename, "c"; attrib=stackmd)
    end
    try
        ncwritevar!(dataset, A)
    finally
        close(dataset)
    end
end

# Stack ########################################################################

struct NCDstack{T,R,W,M,U,K} <: DiskGeoStack{T}
    filename::T
    refdims::R
    window::W
    metadata::M
    usercrs::U
    kwargs::K
end

"""
    NCDstack(filenames; refdims=(), window=(), metadata=nothing)

A lazy GeoStack that loads netcdf files using NCDatasets.jl

Create a stack from a list of filenames.

# Arguments
-`filenames`: `Vector` of `String` paths to netcdf files.

# Keyword arguments
- `refdims`: Add dimension position array was sliced from. Mostly used programatically.
- `window`: can be a tuple of Dimensions, selectors or regular indices.
- `metadata`: Add additional metadata as a `Dict`.
- `keys`: Keys for the layer in each file in filenames. If these do not match a layer
  the first layer will be used. This is also the default.

# Examples
```julia
multifile_stack = NCDstack([path1, path2, path3, path4])
```
"""
NCDstack(filenames::Union{Tuple,Vector};
         refdims=(),
         window=(),
         metadata=nothing,
         keys=cleankeys(ncread(ds -> first(nondimkeys(ds)), fn) for fn in filenames),
         kwargs...) =
    GeoStack(NamedTuple{keys}(filenames), refdims, window, metadata, childtype=NCDarray, kwargs)

"""
    NCDstack(filename; refdims=(), window=(), metadata=nothing)

A lazy GeoStack that loads netcdf files using NCDatasets.jl

Create a stack from the filename of a netcdf file.

# Arguments
-`filename`: `String` path to a netcdf file.

# Keyword arguments
- `refdims`: Add dimension position array was sliced from. Mostly used programatically.
- `window`: can be a tuple of Dimensions, selectors or regular indices.
- `metadata`: Add additional metadata as a `Dict`.

# Examples
```julia
stack = NCDstack(filename; window=(Lat(Between(20, 40),))
stack[:soil_temperature]
```
"""
NCDstack(filename::AbstractString;
         refdims=(),
         window=(),
         metadata=ncread(metadata, filename),
         usercrs=nothing,
         kwargs...) =
    NCDstack(filename, refdims, window, metadata, usercrs, kwargs)

childtype(::NCDstack) = NCDarray
usercrs(stack::NCDarray) = stack.usercrs

# AbstractGeoStack methods

withsource(f, ::Type{NCDarray}, path::AbstractString, key=nothing) = ncread(f, path)
withsourcedata(f, ::Type{NCDarray}, path::AbstractString, key) =
    ncread(d -> f(d[string(key)]), path)

# Override the default to get the dims of the specific key,
# and pass the usercrs from the stack
dims(stack::NCDstack, dataset, key::Key) = dims(dataset, key, usercrs(stack))

missingval(stack::NCDstack) = missing

# Base methods

Base.keys(stack::NCDstack{<:AbstractString}) =
    cleankeys(ncread(nondimkeys, getsource(stack)))

"""
    Base.write(filename::AbstractString, ::Type{NCDstack}, s::AbstractGeoStack)

Write an NCDstack to a single netcdf file, using NCDatasets.jl.

Currently `Dimension` metadata is not handled, and array metadata from other
array types is ignored.
"""
Base.write(filename::AbstractString, ::Type{<:NCDstack}, s::AbstractGeoStack) = begin
    dataset = NCDatasets.Dataset(filename, "c"; attrib=val(metadata(s)))
    try
        map(key -> ncwritevar!(dataset, s[key]), keys(s))
    finally
        close(dataset)
    end
end

# DimensionalData methods for NCDatasets types ###############################

dims(dataset::NCDatasets.Dataset, key::Key, usercrs=nothing) = begin
    v = dataset[string(key)]
    dims = []
    for (i, dimname) in enumerate(NCDatasets.dimnames(v))
        if haskey(dataset, dimname)
            dvar = dataset[dimname]
            # Find the matching dimension constructor. If its an unknown name use
            # the generic Dim with the dim name as type parameter
            dimtype = get(dimmap, dimname, Dim{Symbol(dimname)})
            index = dvar[:]
            mode = _ncdmode(index, dimtype)
            meta = metadata(dvar)

            # Add the dim containing the dimension var array
            push!(dims, dimtype(index, mode, meta))
        else
            # The var doesn't exist. Maybe its `complex` or some other marker,
            # so make it a custom `Dim` with `NoIndex`
            push!(dims, Dim{Symbol(dimname)}(1:size(v, i), NoIndex(), nothing))
        end
    end
    (dims...,)
end

_ncdmode(index, dimtype) = begin
    # Assume the locus is at the center of the cell if boundaries aren't provided.
    # http://cfconventions.org/cf-conventions/cf-conventions.html#cell-boundaries
    # Unless its a time dimension.
    if eltype(index) <: Number
        order = _ncdorder(index)
        span = _ncdspan(index, order)
        sampling = Intervals((dimtype <: TimeDim) ? Start() : Center())
        if dimtype in (Lat, Lon)
            Projected(order, span, sampling, EPSG(4326), nothing)
        else
            Sampled(order, span, sampling)
        end
    elseif eltype(index) <: Dates.AbstractTime
        order = _ncdorder(index)
        span = Irregular(
            if length(index) > 1
                if isrev(indexorder(order))
                    index[end], index[1] + (index[1] - index[2])
                else
                    index[1], index[end] + (index[end] - index[end - 1])
                end
            else
                index[1], index[1]
            end
        )
        sampling = Intervals(Start())
        Sampled(order, span, sampling)
    else
        Categorical()
    end
end

_ncdorder(index) = index[end] > index[1] ? Ordered(Forward(), Forward(), Forward()) :
                                           Ordered(Reverse(), Forward(), Reverse())

_ncdspan(index, order) = begin
    step = index[2] - index[1]
    for i in 2:length(index) -1
        if !(index[i+1] - index[i] ≈ step)
            bounds = if length(index) > 1
                beginhalfcell = abs((index[2] - index[1]) * 0.5)
                endhalfcell = abs((index[end] - index[end-1]) * 0.5)
                if isrev(indexorder(order))
                    index[end] - endhalfcell, index[1] + beginhalfcell
                else
                    index[1] - beginhalfcell, index[end] + endhalfcell
                end
            else
                index[1], index[1]
            end
            return Irregular(bounds)
        end
    end
    return Regular(step)
end


metadata(dataset::NCDatasets.Dataset) = NCDstackMetadata(Dict{String,Any}(dataset.attrib))
metadata(dataset::NCDatasets.Dataset, key::Key) = metadata(dataset[string(key)])
metadata(var::NCDatasets.CFVariable) = NCDarrayMetadata(Dict{String,Any}(var.attrib))
metadata(var::NCDatasets.CFVariable, stackmetadata::NCDstackMetadata) = begin
    md = metadata(var)
    md["dataset"] = stackmetadata
    md
end

missingval(var::NCDatasets.CFVariable) = missing

# https://trac.osgeo.org/gdal/wiki/NetCDF_ProjectionTestingStatus

const cf_proj_params = Dict(
    "false_easting" => "+x_0",
    "false_northing" => "+y_0",
    "scale_factor_at_projection_origin" => "+k_0",
    "scale_factor_at_central_meridian" => "+k_0",
    "standard_parallel[1]" => "+lat_1",
    "standard_parallel[2]" => "+lat_2",
    "longitude_of_central_meridian" => "+lon_0",
    "longitude_of_projection_origin" => "+lon_0",
    "latitude_of_projection_origin" => "+lat_0",
    "straight_vertical_longitude_from_pole" => "+lon_0",
)

const cf_proj_projections = Dict(
    "albers_conical_equal_area" => "+proj=aea",
    "azimuthal_equidistant" => "+proj=aeqd",
    "lambert_azimuthal_equal_area" => "+proj=laea",
    "lambert_conformal_conic" => "+proj=lcc",
    "lambert_cylindrical_equal_area" => "+proj=cea",
    "mercator" => "+proj=merc",
    "orthographic" => "+proj=ortho",
    "polar_stereographic" => "+proj=stere",
    "stereographic" => "+proj=stere",
    "transverse_mercator" => "+proj=tmerc",
)
const proj_cf_projections = Dict(Pair.(values(cf_proj_projections), collect(keys(cf_proj_projections))))


const projections_params = ( 
    albers_conical_equal_area = (
        "standard_parallel[1]",
        "standard_parallel[2]",
        "longitude_of_central_meridian",
        "latitude_of_projection_origin",
        "false_easting",
        "false_northing",
    ),
    azimuthal_equidistant = (
        "longitude_of_projection_origin",
        "latitude_of_projection_origin",
        "false_easting",
        "false_northing",
    ),
    lambert_azimuthal_equal_area = (
        "longitude_of_projection_origin",
        "latitude_of_projection_origin",
        "false_easting",
        "false_northing",
    ),
    lambert_conformal_conic = (
        ("standard_parallel", ["standard_parallel[1]", "standard_parallel[2]"])
        "longitude_of_central_meridian",
        "latitude_of_projection_origin",
        "false_easting",
        "false_northing",
    ),
    lambert_conformal_conic = (
        "longitude_of_central_meridian",
        "latitude_of_projection_origin",
        "false_easting",
        "false_northing",
    ),
    lambert_cylindrical_equal_area = (
        "longitude_of_central_meridian",
        "standard_parallel[1]",
        "false_easting",
        "false_northing",
    ),
    mercator = (
        "longitude_of_projection_origin",
        ("scale_factor_at_projection_origin", "standard_parallel[1]"),
        "false_easting",
        "false_northing",
    ),
    orthographic = (
        "longitude_of_projection_origin",
        "latitude_of_projection_origin",
        "false_easting",
        "false_northing",
    ),
    polar_stereographic = (
        "straight_vertical_longitude_from_pole ",
        "latitude_of_projection_origin",
        ("scale_factor_at_projection_origin", "standard_parallel")
        "false_easting",
        "false_northing",
    ),
    stereographic = (
        "longitude_of_projection_origin",
        "latitude_of_projection_origin",
        "scale_factor_at_projection_origin",
        "false_easting",
        "false_northing",
    ),
    transverse_mercator = (
        "scale_factor_at_central_meridian",
        "longitude_of_central_meridian",
        "latitude_of_projection_origin",
        "false_easting",
        "false_northing",
    )
)

