
# Used internally to expose open disk files inside a `do` block
struct OpenGeoArray{T,N,D<:Tuple,R<:Tuple,A,Na<:Symbol,Me,Mi} <: AbstractGeoArray{T,N,D,A}
    data::A
    dims::D
    refdims::R
    name::Na
    metadata::Me
    missingval::Mi
end
function OpenGeoArray(A, dims, refdims, name, metadata, missingval)
    OpenGeoArray{eltype(A),ndims(A),map(typeof,(dims,refdims,A,name,metadata,missingval))...}(
        A, dims, refdims, name, metadata, missingval
    )
end
function OpenGeoArray(f::Function, A::AbstractGeoArray{T,N}) where {T,N}
    withsourcedata(A) do source
        OA = OpenGeoArray(source, dims(A), refdims(A), name(A), metadata(A), missingval(A))
        f(OA)
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
Base.open(f::Function, A::AbstractGeoArray) = OpenGeoArray(f, A)

@deprecate Open(f, A::AbstractGeoArray) Base.open(f, A)
