# Accept either Symbol or String keys, but allways convert to Symbol
const Key = Union{Symbol,AbstractString}

"""
`AbstractGeoStack` objects hold multipl [`AbstractGeoArray`](@ref) 
that share spatial bounds.

They are `NamedTuple`-like structures that may either contain `NamedTuple`
of [`AbstractGeoArray`](@ref), string paths that will load [`AbstractGeoArray`](@ref), 
or a single path that points to as a file itself containing multiple layers, like 
NetCDF or HDF5. Use and syntax is similar or identical for all cases.

`getindex` on a `AbstractGeoStack` generally returns a memory backed standard 
[`GeoArray`](@ref). `geoarray[:somelayer] |> plot` plots the layers array,
while `geoarray[:somelayer, Lon(1:100), Band(2)] |> plot` will plot the
subset without loading the whole array.

`getindex` on a `AbstractGeoStack` with a key returns another stack with 
getindex applied to all the arrays in the stack.
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


# Interface methods ############################################################

"""
    getsource(s::AbstractGeoStack, [key])

Get the lower-level child object. This can be an `AbstractGeoArray` or
another object with GeoData interface methods defined. 

Working with the low-level object can be better performance as we do not have to
processes everything needed to build a full `AbstractGeoArray`.
"""
function getsource end


"""
    modify(f, series::AbstractGeoStack) 

Apply function `f` to the data of the child `AbstractGeoArray`s. 

`f` must return an idenically sized array.

This method triggers a complete rebuild of all objects, 
and disk based objects will be transferred to memory.
"""
modify(f, stack::AbstractGeoStack) = 
    GeoStack(stack; data=mapdata(A -> modify(f, A), stack))

# Base methods #################################################################

@inline Base.getindex(s::AbstractGeoStack, I...) =
    rebuild(s; data=NamedTuple{keys(s)}(a[I...] for a in values(s)))
# Dict/Array hybrid with dims
@inline Base.getindex(s::AbstractGeoStack, key::Key, I::Vararg{<:Dimension}) =
    getindex(s, key, dims2indices(dims(s, key), I)...)


Base.values(s::AbstractGeoStack) = (s[key] for key in keys(s))
Base.length(s::AbstractGeoStack) = length(keys(s))
Base.keys(s::AbstractGeoStack{<:AbstractString}) = Symbol.(querychild(keys, getsource(s), s))
Base.keys(s::AbstractGeoStack{<:NamedTuple}) = Symbol.(keys(getsource(s)))
Base.names(s::AbstractGeoStack) = keys(s)
Base.map(f, s::AbstractGeoStack) = rebuild(s; data=mapdata(f, s))

mapdata(f, s::AbstractGeoStack) = 
    NamedTuple{cleankeys(keys(s))}(map(f, values(s)))

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
        metadata=metadata(s), childkwargs=childkwargs(s), kwargs...) where T<:MemGeoStack =
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
        metadata=metadata(s), childtype=childtype(s), childkwargs=childkwargs(s), kwargs...) where T<:DiskGeoStack =
    basetypeof(T)(data, refdims, window, metadata, childtype, childkwargs)

getsource(s::DiskGeoStack, args...) = filename(s, args...)

childtype(stack::DiskGeoStack) = stack.childtype

"""
    filename(s::DiskGeoStack)

Return the filename field of a `DiskGeoStack`. 
This may be a `Vector` of `String`, or a `String`.
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
    GeoStack(data...; keys, kwargs...)
    GeoStack(data::Union{Vector,Tuple}; keys, kwargs...)
    GeoStack(data::NamedTuple; window=(), metadata=nothing, refdims=(), childkwargs=()) =
    GeoStack(s::AbstractGeoStack; [keys, data, refdims, window, metadata])

A concrete `MemGeoStack` implementation. Holds layers of [`GeoArray`](@ref).

## Argumenst

- `data`: A `NamedTuple` of [`GeoArray`](@ref), or a `Vector`, `Tuple` or splatted arguments
  of [`GeoArray`](@ref). The latter options must pass a `keys` keyword argument.

## Keyword Argumenst

- `keys`: Used as stack keys when a `Tuple` or `Vector` or splat of geoarrays are passed in.
- `window`: A `Tuple` of `Dimension`/`Selector`/indices that will be applied to the 
  contained arrays when they are accessed.
- `refdims`: Reference dimensions from earlier subsetting.
- `metadata`: Metadata as a [`StackMetadata`](@ref) object.
- `childkwargs`: A `NamedTuple` of keyword arguments to pass to the constructor.
- `refdims`: `Tuple` of  position `Dimension` the array was sliced from.
"""
struct GeoStack{T,R,W,M,K} <: MemGeoStack{T}
    data::T
    refdims::R
    window::W
    metadata::M
    childkwargs::K
end
GeoStack(data::AbstractGeoArray...; keys=name.(data), kwargs...) =
    GeoStack(NamedTuple{cleankeys(keys)}(data); kwargs...)
GeoStack(data::NamedTuple; 
         refdims=(), 
         window=(), 
         metadata=nothing, 
         childkwargs=()) =
    GeoStack(data, refdims, window, metadata, childkwargs)
GeoStack(s::AbstractGeoStack; keys=cleankeys(Base.keys(s)),
         data=NamedTuple{keys}(s[key] for key in keys),
         refdims=refdims(s), 
         window=(), 
         metadata=metadata(s), 
         childkwargs=()) =
    GeoStack(data, refdims, window, metadata, childkwargs)

Base.convert(::Type{GeoStack}, src::AbstractGeoStack) = GeoStack(src)



# Concrete DiskGeoStack implementation ######################################################

"""
    DiskStack(filenames...; keys, kwargs...)
    DiskStack(filenames; keys, kwargs...)
    DiskStack(filenames::NamedTuple; 
              window=(), 
              metadata=nothing, 
              childtype, 
              childkwargs=()
              refdims=())

Concrete [`DiskGeoStack`](@ref) implementation. Loads a stack of files lazily from disk.

## Arguments

- `filename`: a NamedTuple of stack keys and `String` filenames.

## Keyword arguments

- `keys`: Used as stack keys when a `Tuple`, `Vector` or splat of filenames are passed in.
- `window`: A `Tuple` of `Dimension`/`Selector`/indices that will be applied to the 
  contained arrays when they are accessed.
- `metadata`: Metadata as a [`StackMetadata`](@ref) object.
- `childtype`: The type of the child data. eg. `GDALarray`. Required.
- `childkwargs`: A `NamedTuple` of keyword arguments to pass to the `childtype` constructor.
- `refdims`: `Tuple` of  position `Dimension` the array was sliced from.
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
DiskStack(filenames...; kwargs...) = DiskStack(filenames; kwargs...)
DiskStack(filenames; keys, kwargs...) =
    DiskStack(NamedTuple{cleankeys(keys)}((filenames...,)); kwargs...)


# Other base methods

Base.copy(stack::AbstractGeoStack) =
    rebuild(stack; data=map(copy, getsource(stack)))

"""
    Base.copy!(dst::AbstractGeoStack, src::AbstractGeoStack, [keys=keys(dst)])

Copy all or a subset of layers from one stack to another.

## Example

Copy just the `:sea_surface_temp` and `:humidity` layers from `src` to `dst`.

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

## Example

Copy the `:humidity` layer from `stack` to `array`.

```julia
copy!(array, stack, :humidity)
```
"""
Base.copy!(dst::AbstractArray, src::AbstractGeoStack, key) =
    copy!(dst, src[key])

"""
    Base.cat(stacks::AbstractGeoStack...; [keys=keys(stacks[1])], dims)

Concatenate all or a subset of layers for all passed in stacks.

## Keyword Arguments

- `keys`: `Tuple` of `Symbol` for the stack keys to concatenate.
- `dims`: Dimension of child array to concatenate on.

## Example

Concatenate the :sea_surface_temp and :humidity layers in the time dimension:

```julia
cat(stacks...; keys=(:sea_surface_temp, :humidity), dims=Ti)
```
"""
Base.cat(stacks::AbstractGeoStack...; keys=keys(stacks[1]), dims) = begin
    vals = Tuple(cat((s[key] for s in stacks)...; dims=dims) for key in keys)
    GeoStack(stacks[1], data=NamedTuple{keys}(vals))
end
Base.first(s::AbstractGeoStack) = s[first(keys(s))]
Base.last(s::AbstractGeoStack) = s[last(keys(s))]
