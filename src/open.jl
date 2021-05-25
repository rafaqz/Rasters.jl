
"""
    OpenGeoArray <: AbstractGeoArray

Used internally to expose open disk files inside a `do` block
"""
struct OpenGeoArray{T,N,D<:Tuple,R<:Tuple,A<:AbstractArray{T,N},Na<:Symbol,Me,Mi} <: AbstractGeoArray{T,N,D,A}
    data::A
    dims::D
    refdims::R
    name::Na
    metadata::Me
    missingval::Mi
end
function OpenGeoArray(f::Function, A::AbstractGeoArray{T,N}; kw...) where {T,N}
    # Open FileArray to expose the actual dataset object, even inside nested wrappers
    fa = Flatten.flatten(data(A), FileArray)
    if fa == ()
        f(OpenGeoArray(data(A), dims(A), refdims(A), name(A), metadata(A), missingval(A)))
    else
        open(fa[1]; kw...) do x
            # Rewrap the opened object where the FileArray was
            d = Flatten.reconstruct(data(A), (x,), FileArray) 
            f(OpenGeoArray(d, dims(A), refdims(A), name(A), metadata(A), missingval(A)))
        end
    end
end

"""
    open(f, A::AbstractGeoArray)

`open` is used to open any `AbstractGeoArray` and do multiple operations
on it in a safe way. It's a shorthand for the unexported `OpenGeoArray`
constructor.

`f` is a method that accepts a single argument - an `OpenGeoArray` object
which is just an `AbstractGeoArray` that holds an open disk - based object.
Often it will be a `do` block:

```julia
ga = GDALarray(filepath)
open(ga) do A
    A[I...] # A is an `OpenGeoArray` wrapping the disk-based object.
    # ...  multiple things you need to do with the open file
end
```

By using a do block to open file we ensure they are always closed again
after we finish working with them.
"""
Base.open(f::Function, A::AbstractGeoArray; write=false) = OpenGeoArray(f, A; write)

@deprecate Open(f, A::AbstractGeoArray) Base.open(f, A)
