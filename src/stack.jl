# Accept either Symbol or String keys, but allways convert to Symbol
const Key = Union{Symbol,AbstractString}

"""
Stack objects hold multiple raster array that share spatial metadata and bounds.

These are `NamedTuple`-like structures that may either contain `NamedTuple`
of `AbstractGeoArray`, string paths that will load `AbstractGeoArray`, or a single
path that points to as a file structured as a multi-layered stack, like NetCDF.

The primary purpose of  is that use and syntax is identical for all cases,
abstracting away data source and simplifying access code. `getindex` on any
`AbstractGeoStack` may return a memory backed standard `GeoArray`, or a disk
base AbstractGeoArray. `geoarray[:somelayer] |> plot` plots the layers array,
while `geoarray[:somelayer, Lon(1:100), Band(2)] |> plot` will plot the
subsetted array directly from disk, without loading the whole array.

`getindex` or `view` on a GeoStack returns another stack with the method applied
to all the arrays in the stack.
"""
abstract type AbstractGeoStack{T} end

(::Type{T})(data, keys; kwargs...) where T<:AbstractGeoStack =
    T(NamedTuple{Tuple(Symbol.(keys))}(Tuple(data)); kwargs...)

# Standard fields

refdims(s::AbstractGeoStack) = s.refdims
metadata(s::AbstractGeoStack) = s.metadata
window(s::AbstractGeoStack) = s.window

"""
    childkwargs(s::AbstractGeoStack)

Returns the keyword arguments that will be passed to the child array constructor.
"""
childkwargs(s::AbstractGeoStack) = s.childkwargs

@inline Base.getindex(s::AbstractGeoStack, I...) =
    rebuild(s; data=NamedTuple{keys(s)}(a[I...] for a in values(s)))

# Interface methods ############################################################


"""
    getsource(s::AbstractGeoStack, [key])

Get the lower lovel child object. This can be an `AbstractGeoArray` or
a lower-level object with GeoData methods defined. Returning
the low-level object can be better performance as we do not have to
processes everything needed to build a full `AbstractGeoArray`.
"""
function getsource end

# Base methods #################################################################

# Dict/Array hybrid with dims
@inline Base.getindex(s::AbstractGeoStack, key::Key, I::Vararg{<:Dimension}) =
    getindex(s, key, dims2indices(dims(s, key), I)...)

Base.values(s::AbstractGeoStack) = (s[key] for key in keys(s))
Base.length(s::AbstractGeoStack) = length(keys(s))
Base.keys(s::AbstractGeoStack{<:AbstractString}) = Symbol.(querychild(keys, getsource(s), s))
Base.keys(s::AbstractGeoStack{<:NamedTuple}) = Symbol.(keys(getsource(s)))
Base.names(s::AbstractGeoStack) = keys(s)

"""
    Base.write(filename::AbstractString, T::Type{<:AbstractGeoArray}, s::AbstractGeoStack)

Save all layers of an `AbstractGeoStack` to separate files, using the backend determined
by `T`.

## Example

```julia
write(filename, GDALarray, A)
```
"""
Base.write(filename::AbstractString, ::Type{T}, s::AbstractGeoStack) where T <: AbstractGeoArray =
    for key in keys(s)
        base, ext = splitext(filename)
        fn = joinpath(string(base, "_", key, ext))
        write(fn, T, s[key])
    end


# Memory-based stacks ######################################################

"""
An [`AbstractGeoStack`](@ref) stored in memory.
"""
abstract type MemGeoStack{T} <: AbstractGeoStack{T} end

data(s::MemGeoStack) = s.data
data(s::MemGeoStack{<:NamedTuple}, key::Key) = data(s)[Symbol(key)]

rebuild(s::T; data=data(s), refdims=refdims(s), window=window(s),
        metadata=metadata(s), childkwargs=childkwargs(s)) where T<:MemGeoStack =
    basetypeof(T)(data, refdims, window, metadata, childkwargs)

getsource(s::MemGeoStack{<:NamedTuple}, args...) = data(s, args...)

childdata(f, childobj, ::AbstractGeoStack) = f(childobj)

@inline Base.view(s::MemGeoStack, I...) =
    rebuild(s; data=NamedTuple{keys(s)}(view(a, I...) for a in values(s)))

@inline Base.getindex(s::MemGeoStack, key::Key) = begin
    A = rebuild(data(s)[key]; childkwargs(s)...)
    window_ = maybewindow2indices(dims(A), window(s))
    readwindowed(A, window_)
end
@inline Base.getindex(s::MemGeoStack, key::Key, i1::StandardIndices, I::StandardIndices...) = begin
    A = rebuild(data(s)[key]; childkwargs(s)...)
    window_ = maybewindow2indices(dims(A), window(s))
    readwindowed(A, window_, i1, I...)
end

# Disk-based stacks ######################################################

"""
[`AbstractGeoStack`](@ref)s stored on disk.
"""
abstract type DiskGeoStack{T} <: AbstractGeoStack{T} end

rebuild(s::T; data=filename(s), refdims=refdims(s), window=window(s),
        metadata=metadata(s), childtype=childtype(s), childkwargs=childkwargs(s)) where T<:DiskGeoStack =
    basetypeof(T)(data, refdims, window, metadata, childtype, childkwargs)

getsource(s::DiskGeoStack, args...) = filename(s, args...)

childtype(stack::DiskGeoStack) = stack.childtype

"""
    filename(s::DiskGeoStack)

Return the filename field of a `DiskGeoStack`. This may be a `Vector` of `String`,
or a `String`.
"""
filename(s::DiskGeoStack) = s.filename
"""
    filename(s::DiskGeoStack, key)

Return the filename field of a `DiskGeoStack` for a given key.

This will always be a single string. However, in some cases
all keys may have the same filename.
"""
filename(s::DiskGeoStack{<:NamedTuple}, key::Key) = filename(s)[Symbol(key)]
filename(s::DiskGeoStack, key::Key) = filename(s)


# Default to using the same method as GeoStack.
# But a method for withsourcedata dispatching on the
# array object may save time in some cases, like when
# the array size has to be determined before loading the data
# like with mmaped grd files.
withsourcedata(f, A::DiskGeoArray, key...) =
    withsourcedata(f, typeof(A), filename(A), key...)
# By default assume the source object is a single object with
# all info and data.
withsourcedata(f, childtype::Type, path, key...) =
    withsource(f, childtype, path, key...)

withsource(f, childtype::Type, path, key...) = f(path)

# Base methods

Base.getindex(s::DiskGeoStack, key::Key) = begin
    filename_ = filename(s, key)
    withsource(childtype(s), filename_, key) do dataset
        A = childtype(s)(dataset, filename_, key; name=string(key), childkwargs(s)...)
        window_ = maybewindow2indices(A, window(s))
        readwindowed(A, window_)
    end
end
Base.getindex(s::DiskGeoStack, key::Key, i1::StandardIndices, I::StandardIndices...) = begin
    filename_ = filename(s, key)
    withsource(childtype(s), filename_, key) do dataset
        A = childtype(s)(dataset, filename_, key; name=string(key), childkwargs(s)...)
        window_ = maybewindow2indices(dims(A), window(s))
        readwindowed(A, window_, i1, I...)
    end
end

@inline Base.view(s::DiskGeoStack, I...) = rebuild(s; window=I)


# Default dims, metadata and missingval methods
#
# For DiskStack we query the underlying object - avoiding building
# an AbstractGeoArray unless we have to. Examples are be an NCDatasets,
# ArchGDAL, or HDF5 `Dataset` object. These sources will add a
# For MemGeoStack the childobj is just a GeoArray as it is allready loaded in memory.
for func in (:dims, :metadata, :missingval)
    @eval begin
        $func(s::AbstractGeoStack, key::Key) = $func(s[key])
        $func(s::DiskGeoStack, key::Key) =
            withsource(childtype(s), getsource(s, key), key) do sourceobj
                $func(sourceobj, key)
            end
    end
end

# Concrete MemGeoStack implementation ######################################################

"""
Concrete `MemGeoStack` implementation. Holds concrete [`GeoArray`](@ref) layers in memory.

- `data`: A `NamedTuple` of `GeoArray`.
- `window`: A `Tuple` of Dimensions/Selectors/Indices that will be applied to the contained
            arrays when they are accessed.
- `refdims`: reference dimensions from earlier subsetting.
- `metadata`: Any metadata object.
"""
struct GeoStack{T,R,W,M,K} <: MemGeoStack{T}
    data::T
    refdims::R
    window::W
    metadata::M
    childkwargs::K
end
"""
    GeoStack(data::Vararg{<:AbstractGeoArray}; keys, childkwargs)

Convert `GeoArray`s to a `GeoStack`.
"""
GeoStack(data::AbstractGeoArray...; keys=name.(data), childkwargs=()) =
    GeoStack(NamedTuple{cleankeys(keys)}(data); childkwargs=())
"""
    GeoStack(data::NamedTuple; [window=()], [metadata=nothing], childkwargs) =

Construct a `GeoStack` from a NamedTuple of [`GeoArray`](@ref) and keyword arguments.

The `childkwargs` keyword is used as keyword arguments for the child contructor.
"""
GeoStack(data::NamedTuple; refdims=(), window=(), metadata=nothing, childkwargs=()) =
    GeoStack(data, refdims, window, metadata, childkwargs)
"""
    GeoStack(s::AbstractGeoStack; [keys, data, refdims, window, metadata])

Construct a `GeoStack` from another `AbstractGeoStack` and keyword arguments.
`data` must be a `NamedTuple` of `GeoArray`.
"""
GeoStack(s::AbstractGeoStack; keys=cleankeys(Base.keys(s)),
         data=NamedTuple{keys}(s[key] for key in keys),
         refdims=refdims(s), window=(), metadata=metadata(s), childkwargs=()) =
    GeoStack(data, refdims, window, metadata, childkwargs)

Base.convert(::Type{GeoStack}, src::AbstractGeoStack) = GeoStack(src)



# Concrete DiskGeoStack implementation ######################################################

"""
    DiskStack(filename::NamedTuple; window=())

Load a stack of files lazily from disk.

# Arguments
- `filename`: a NamedTuple of `String` filenames.

# Keyword arguments
- `window`: can be a tuple of Dimensions, selectors or regular indices.
- `childtype`: the type of the child data. eg. `GDALarray`.
- `kwargs`: keyword arguments to pass to the child constructor
"""
struct DiskStack{T,R,W,M,C,K} <: DiskGeoStack{T}
    filename::T
    refdims::R
    window::W
    metadata::M
    childtype::C
    childkwargs::K
end
DiskStack(filenames::NamedTuple; refdims=(), window=(), metadata=nothing, childtype, childkwargs=()) =
    DiskStack(filenames, refdims, window, metadata, childtype, childkwargs)
DiskStack(filenames, keys; kwargs...) =
    DiskStack(NamedTuple{cleankeys(keys)}((filenames...,)); kwargs...)


# Other base methods

Base.copy(stack::AbstractGeoStack) =
    rebuild(stack; data=map(copy, getsource(stack)))

"""
    Base.copy!(dst::AbstractGeoStack, src::AbstractGeoStack, [keys=keys(dst)])

Copy all or a subset of layers from one stack to another.

## Example
```julia
copy!(dst::AbstractGeoStack, src::AbstractGeoStack, keys=(:sea_surface_temp, :humidity))
```
"""
Base.copy!(dst::MemGeoStack, src::AbstractGeoStack, keys=keys(dst)) = begin
    for key in keys
        key in Symbol.(Base.keys(dst)) || throw(ArgumentError("key $key not found in dest keys"))
        key in Symbol.(Base.keys(src)) || throw(ArgumentError("key $key not found in source keys"))
    end
    for key in Symbol.(keys)
        copy!(dst[key], src[key])
    end
end
"""
    Base.copy!(dst::AbstractArray, src::DiskGeoStack, key::Key)

Copy the stack layer `key` to `dst`, which can be any `AbstractArray`.
"""
Base.copy!(dst::AbstractArray, src::AbstractGeoStack, key) =
    copy!(dst, src[key])

"""
    Base.cat(stacks::AbstractGeoStack...; [keys=keys(stacks[1])], dims)

Concatenate all or a subset of layers for all passed in stacks.

## Example
```julia
cat(stacks...; keys=(:sea_surface_temp, :humidity), dims)
```
"""
Base.cat(stacks::AbstractGeoStack...; keys=keys(stacks[1]), dims) = begin
    vals = Tuple(cat((s[key] for s in stacks)...; dims=dims) for key in keys)
    GeoStack(stacks[1], data=NamedTuple{keys}(vals))
end
Base.first(s::AbstractGeoStack) = s[first(keys(s))]
Base.last(s::AbstractGeoStack) = s[last(keys(s))]
